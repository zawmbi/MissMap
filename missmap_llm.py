"""
Provider-agnostic LLM harness for synonym recovery (Aim 3).

Design rules that make the benchmark defensible:
  * one versioned prompt (PROMPT_VERSION) for every model;
  * every call is logged to JSONL with the raw response, exact model ID,
    mode, run index, latency and token usage;
  * a call with the same (provider, model, mode, prompt version, species,
    run index) is served from the log instead of re-billed;
  * failures are recorded as failures (error field), never as "no synonyms";
  * scientific names and common names come back in separate fields, so
    common names can't inflate synonym counts.

Model spec strings (see make_provider):
    anthropic:claude-opus-5-5           closed-book (memory only)
    anthropic:claude-sonnet-5-5+web     with Anthropic's server-side web search
    groq:llama-3.3-70b-versatile        Groq (OpenAI-compatible endpoint)
    openai:<model>                      OpenAI
    openrouter:<vendor/model>           OpenRouter (Qwen, DeepSeek, Gemini, Mistral ...)
    ollama:<model>                      local Ollama server (no cost, fully reproducible weights)
"""

import json
import os
import re
import time
from dataclasses import dataclass

import requests

PROMPT_VERSION = "v2-2026-11-19"

SYSTEM_PROMPT = (
    "You are an expert in botanical and zoological nomenclature. "
    "Given one accepted species name, list every scientific name that a major "
    "taxonomic authority treats, or has treated, as a synonym of that species: "
    "basionyms, homotypic and heterotypic synonyms, and orthographic variants. "
    "Give each name without authorship in `name` and with the author citation in "
    "`authorship`. Do not include the input name itself, cultivars, hybrids, "
    "infraspecific names, other accepted species, or common names in `synonyms`; "
    "put vernacular names only in `common_names`. If you are not confident a name "
    "exists, leave it out. If you know of none, return an empty list."
)

USER_TEMPLATE = (
    "Species: {species}\n"
    "{known_block}"
    "Return JSON matching the schema."
)

RESPONSE_SCHEMA = {
    "type": "object",
    "properties": {
        "accepted_name": {"type": "string"},
        "synonyms": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "authorship": {"type": "string"},
                    "type": {"type": "string", "enum": [
                        "basionym", "homotypic", "heterotypic",
                        "orthographic_variant", "misapplied", "unknown"]},
                    "source": {"type": "string"},
                },
                "required": ["name", "authorship", "type", "source"],
                "additionalProperties": False,
            },
        },
        "common_names": {"type": "array", "items": {"type": "string"}},
    },
    "required": ["accepted_name", "synonyms", "common_names"],
    "additionalProperties": False,
}

WEB_SEARCH_SUFFIX = (
    "\nYou may search the web. Prefer POWO/WCVP, World Flora Online, IPNI, "
    "Tropicos, GBIF, Catalogue of Life and recent taxonomic revisions. "
    "Finish with ONLY the JSON object (no prose) matching this schema:\n"
    + json.dumps(RESPONSE_SCHEMA)
)


def build_user_message(species, known=None):
    known_block = ""
    if known:
        known_block = ("Already known (do not repeat): " + json.dumps(sorted(set(known))) + "\n")
    return USER_TEMPLATE.format(species=species, known_block=known_block)


def extract_json(text):
    """Parse the JSON object in a model reply (tolerates code fences / prose around it)."""
    if not text:
        raise ValueError("empty response")
    t = text.strip()
    m = re.search(r"```(?:json)?\s*(\{.*\})\s*```", t, re.DOTALL)
    if m:
        t = m.group(1)
    else:
        start, end = t.find("{"), t.rfind("}")
        if start == -1 or end == -1:
            raise ValueError("no JSON object in response")
        t = t[start:end + 1]
    return json.loads(t)


def normalise_result(parsed):
    """Coerce a parsed reply into the schema; accepts the legacy {'synonyms': [str]} shape."""
    if not isinstance(parsed, dict):
        raise ValueError("response is not a JSON object")
    syns = []
    for s in parsed.get("synonyms", []) or []:
        if isinstance(s, str):
            syns.append({"name": s.strip(), "authorship": "", "type": "unknown", "source": ""})
        elif isinstance(s, dict) and str(s.get("name", "")).strip():
            syns.append({"name": str(s["name"]).strip(),
                         "authorship": str(s.get("authorship", "") or "").strip(),
                         "type": str(s.get("type", "unknown") or "unknown"),
                         "source": str(s.get("source", "") or "")})
    return {"accepted_name": str(parsed.get("accepted_name", "") or ""),
            "synonyms": syns,
            "common_names": [str(c) for c in parsed.get("common_names", []) or [] if c]}


@dataclass
class CallResult:
    provider: str
    model: str
    mode: str
    species: str
    run_index: int
    prompt_version: str
    timestamp: str
    latency_s: float
    raw_text: str
    parsed: dict
    error: str
    usage: dict
    settings: dict

    def to_json(self):
        return json.dumps(self.__dict__, ensure_ascii=False)


# ----------------------------------------------------------------------
# Providers
# ----------------------------------------------------------------------
class AnthropicProvider:
    """
    Claude via the official Anthropic SDK (needs an API key or `ant auth login`;
    a claude.ai Pro/Max subscription does not cover API calls).

    Notes for the benchmark:
      * Current Claude models reject `temperature`; run-to-run variation is
        measured with repeated runs instead.
      * Thinking cannot be disabled on Claude Opus 5.5; depth is set with
        `effort`, which is recorded in `settings`.
      * Server-side refusal fallbacks are deliberately NOT enabled: a fallback
        would silently answer with a different model and contaminate the
        per-model scores. A refusal is logged as an error.
    """

    name = "anthropic"

    def __init__(self, model, web_search=False, effort=None, max_web_searches=5, client=None):
        self.model = model
        self.web_search = web_search
        self.effort = effort
        self.max_web_searches = max_web_searches
        self._client = client

    @property
    def mode(self):
        return "web" if self.web_search else "closed"

    @property
    def settings(self):
        return {"effort": self.effort, "web_search": self.web_search,
                "max_web_searches": self.max_web_searches if self.web_search else None}

    def client(self):
        if self._client is None:
            import anthropic  # imported lazily so the rest of MissMap runs without it
            self._client = anthropic.Anthropic()
        return self._client

    def complete(self, species, known=None):
        kwargs = {"model": self.model, "max_tokens": 16000,
                  "system": SYSTEM_PROMPT + (WEB_SEARCH_SUFFIX if self.web_search else "")}
        output_config = {}
        if self.effort:
            output_config["effort"] = self.effort
        messages = [{"role": "user", "content": build_user_message(species, known)}]
        if self.web_search:
            kwargs["tools"] = [{"type": "web_search_20260209", "name": "web_search",
                                "max_uses": self.max_web_searches}]
        else:
            # Structured outputs guarantee schema-valid JSON in closed-book mode.
            # (Not combined with web search: search results carry citations,
            # which structured outputs do not accept.)
            output_config["format"] = {"type": "json_schema", "schema": RESPONSE_SCHEMA}
        if output_config:
            kwargs["output_config"] = output_config

        usage = {"input_tokens": 0, "output_tokens": 0, "web_search_requests": 0}
        text_parts = []
        for _ in range(5):  # resume server-tool turns that pause
            resp = self.client().messages.create(messages=messages, **kwargs)
            u = resp.usage
            usage["input_tokens"] += getattr(u, "input_tokens", 0) or 0
            usage["output_tokens"] += getattr(u, "output_tokens", 0) or 0
            stu = getattr(u, "server_tool_use", None)
            if stu is not None:
                usage["web_search_requests"] += getattr(stu, "web_search_requests", 0) or 0
            if resp.stop_reason == "refusal":
                cat = getattr(getattr(resp, "stop_details", None), "category", None)
                raise RuntimeError(f"refusal (category={cat})")
            text_parts = [b.text for b in resp.content if b.type == "text"]
            if resp.stop_reason == "pause_turn":
                messages = messages + [{"role": "assistant", "content": resp.content}]
                continue
            if resp.stop_reason == "max_tokens":
                raise RuntimeError("hit max_tokens before finishing")
            break
        usage["request_id"] = getattr(resp, "_request_id", None)
        # With web search, citations split the reply into several text blocks; join them.
        return "".join(text_parts), usage


class OpenAICompatProvider:
    """Any OpenAI-style /chat/completions endpoint (Groq, OpenAI, OpenRouter, Together, Ollama)."""

    ENDPOINTS = {
        "groq": ("https://api.groq.com/openai/v1", "GROQ_API_KEY"),
        "openai": ("https://api.openai.com/v1", "OPENAI_API_KEY"),
        "openrouter": ("https://openrouter.ai/api/v1", "OPENROUTER_API_KEY"),
        "together": ("https://api.together.xyz/v1", "TOGETHER_API_KEY"),
        "ollama": ("http://localhost:11434/v1", None),
    }

    def __init__(self, provider, model, temperature=0.0, json_mode=True, base_url=None, api_key=None):
        if provider not in self.ENDPOINTS and not base_url:
            raise ValueError(f"Unknown provider '{provider}'")
        self.name = provider
        self.model = model
        self.temperature = temperature
        self.json_mode = json_mode
        default_url, key_env = self.ENDPOINTS.get(provider, (base_url, None))
        self.base_url = base_url or default_url
        self.api_key = api_key or (os.getenv(key_env) if key_env else None)
        if key_env and not self.api_key:
            raise KeyError(f"Missing {key_env} for provider '{provider}'")

    mode = "closed"

    @property
    def settings(self):
        return {"temperature": self.temperature, "json_mode": self.json_mode}

    def complete(self, species, known=None):
        payload = {
            "model": self.model,
            "temperature": self.temperature,
            "max_tokens": 4000,
            "messages": [
                {"role": "system", "content": SYSTEM_PROMPT + "\nSchema: " + json.dumps(RESPONSE_SCHEMA)},
                {"role": "user", "content": build_user_message(species, known)},
            ],
        }
        if self.json_mode:
            payload["response_format"] = {"type": "json_object"}
        headers = {"Content-Type": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        last = None
        for attempt in range(4):
            r = requests.post(f"{self.base_url}/chat/completions", headers=headers,
                              json=payload, timeout=120)
            if r.status_code in (429, 500, 502, 503, 504):
                last = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
                time.sleep(2.0 * 2 ** attempt)
                continue
            if r.status_code != 200:
                raise RuntimeError(f"HTTP {r.status_code}: {r.text[:300]}")
            data = r.json()
            u = data.get("usage", {}) or {}
            return (data["choices"][0]["message"]["content"] or "",
                    {"input_tokens": u.get("prompt_tokens"), "output_tokens": u.get("completion_tokens"),
                     "served_model": data.get("model")})
        raise last


def make_provider(spec, temperature=0.0, effort=None):
    """Build a provider from 'provider:model[+web]'."""
    if ":" not in spec:
        raise ValueError(f"Model spec must look like 'provider:model', got '{spec}'")
    provider, model = spec.split(":", 1)
    web = model.endswith("+web")
    model = model[:-4] if web else model
    if provider == "anthropic":
        return AnthropicProvider(model, web_search=web, effort=effort)
    if web:
        raise ValueError("'+web' is only wired up for anthropic models")
    return OpenAICompatProvider(provider, model, temperature=temperature)


# ----------------------------------------------------------------------
# Logged, cached calls
# ----------------------------------------------------------------------
class LLMLog:
    """Append-only JSONL log that doubles as a cache of successful calls."""

    def __init__(self, path):
        self.path = path
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self._done = {}
        if os.path.exists(path):
            with open(path, encoding="utf-8") as f:
                for line in f:
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if not rec.get("error"):
                        self._done[self.key(rec["provider"], rec["model"], rec["mode"],
                                            rec["prompt_version"], rec["species"], rec["run_index"])] = rec

    @staticmethod
    def key(provider, model, mode, prompt_version, species, run_index):
        return f"{provider}|{model}|{mode}|{prompt_version}|{species.lower()}|{run_index}"

    def get(self, provider, species, run_index):
        return self._done.get(self.key(provider.name, provider.model, provider.mode,
                                       PROMPT_VERSION, species, run_index))

    def append(self, result):
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(result.to_json() + "\n")
        if not result.error:
            self._done[self.key(result.provider, result.model, result.mode, result.prompt_version,
                                result.species, result.run_index)] = json.loads(result.to_json())

    def records(self):
        if not os.path.exists(self.path):
            return []
        with open(self.path, encoding="utf-8") as f:
            return [json.loads(l) for l in f if l.strip()]


def query_synonyms(provider, species, known=None, run_index=0, log=None):
    """One logged call. Returns a dict (the JSONL record)."""
    if log is not None:
        cached = log.get(provider, species, run_index)
        if cached is not None:
            return cached
    t0 = time.time()
    raw, usage, parsed, err = "", {}, None, ""
    try:
        raw, usage = provider.complete(species, known)
        parsed = normalise_result(extract_json(raw))
    except Exception as e:  # noqa: BLE001 - recorded, not swallowed
        err = f"{type(e).__name__}: {e}"
    res = CallResult(provider=provider.name, model=provider.model, mode=provider.mode,
                     species=species, run_index=run_index, prompt_version=PROMPT_VERSION,
                     timestamp=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                     latency_s=round(time.time() - t0, 3), raw_text=raw, parsed=parsed,
                     error=err, usage=usage, settings=provider.settings)
    if log is not None:
        log.append(res)
    return json.loads(res.to_json())


def synonym_names(record):
    """Scientific synonym names from a logged record ([] when the call failed)."""
    p = record.get("parsed") or {}
    return [s["name"] for s in p.get("synonyms", []) if s.get("name")]
