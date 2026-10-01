"""The one door to the network.

Every HTTP request the package makes goes through :func:`urlopen`, under a
:class:`NetworkPolicy` taken from ``knowledge.network``:

- ``offline`` (the default) refuses every request;
- ``institution`` allows localhost and the hosts in ``knowledge.allowed_hosts``,
  so a local or institutional LLM, ontology mirror, or table server can be used
  while nothing leaves the institution;
- ``online`` allows any host.

Redirects are re-checked, so an allowed host cannot bounce a request outside.
Census uses its own client (tiledb), so plan offers a Census route only in
``online`` mode. CellTypist routes read cached models and never download.
"""

from __future__ import annotations

import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any, Literal

from .contracts import ContractError

NetworkMode = Literal["offline", "institution", "online"]
LOOPBACK = frozenset({"localhost", "127.0.0.1", "::1"})


@dataclass(frozen=True)
class NetworkPolicy:
    mode: NetworkMode = "offline"
    allowed_hosts: tuple[str, ...] = ()

    def check(self, url: str) -> None:
        parsed = urllib.parse.urlsplit(url)
        host = (parsed.hostname or "").lower()
        if parsed.scheme not in {"http", "https"}:
            raise ContractError(f"only http(s) requests are allowed, not {url!r}")
        if self.mode == "online":
            return
        if self.mode == "offline":
            raise ContractError(
                f"network is offline (knowledge.network); refusing {host or url!r}"
            )
        # An entry starting with "." allows that domain's subdomains.
        allowed = host in LOOPBACK or any(
            host == entry or (entry.startswith(".") and host.endswith(entry))
            for entry in (item.lower() for item in self.allowed_hosts)
        )
        if not allowed:
            raise ContractError(
                f"{host!r} is outside the institution (knowledge.allowed_hosts); "
                "set knowledge.network: online to allow it"
            )


OFFLINE = NetworkPolicy()
# Only for explicit download commands (acquire-*) that a person runs by name and
# that send nothing but a GET for a public, checksum-pinned file.
ONLINE = NetworkPolicy("online")


class _CheckedRedirect(urllib.request.HTTPRedirectHandler):
    def __init__(self, policy: NetworkPolicy, follow: bool) -> None:
        self.policy = policy
        self.follow = follow

    def redirect_request(
        self, req: Any, fp: Any, code: int, msg: str, headers: Any, newurl: str
    ) -> Any:
        # urllib copies headers, Authorization included, to the redirect target.
        if not self.follow:
            raise ContractError(f"refusing a {code} redirect to {newurl!r}")
        self.policy.check(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open(
    request: urllib.request.Request, timeout: float, policy: NetworkPolicy, follow: bool = True
) -> Any:
    opener = urllib.request.build_opener(_CheckedRedirect(policy, follow))
    return opener.open(request, timeout=timeout)


def urlopen(
    request: urllib.request.Request,
    *,
    timeout: float,
    policy: NetworkPolicy,
    redirects: bool = True,
) -> Any:
    policy.check(request.full_url)
    return _open(request, timeout, policy, redirects)
