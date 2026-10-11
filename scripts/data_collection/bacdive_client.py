# Shared BacDive API helpers (base URLs, JSON fetch with one User-Agent and timeout) for the curation scripts.

from __future__ import annotations

import json
from typing import Any
from urllib.request import Request, urlopen


BACDIVE_API_BASE = "https://api.bacdive.dsmz.de/v2"
BACDIVE_PORTAL_BASE = "https://bacdive.dsmz.de/strain"
BACDIVE_USER_AGENT = "extremophile-thesis-curation/1.0"
BACDIVE_TIMEOUT_SECONDS = 30


def ensure_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def portal_url(bacdive_id: int | str) -> str:
    return f"{BACDIVE_PORTAL_BASE}/{bacdive_id}"


def fetch_json(
    url: str,
    *,
    user_agent: str = BACDIVE_USER_AGENT,
    timeout: int = BACDIVE_TIMEOUT_SECONDS,
) -> dict[str, Any]:
    # All BacDive HTTP requests flow through one place so headers/timeouts stay
    # consistent across the different curation scripts.
    request = Request(url, headers={"User-Agent": user_agent})
    with urlopen(request, timeout=timeout) as response:
        return json.load(response)
