from __future__ import annotations

import urllib.request

import pytest

from cell_curator.config import KnowledgeConfig
from cell_curator.contracts import ContractError
from cell_curator.net import NetworkPolicy, _CheckedRedirect, urlopen


def test_each_mode_reaches_exactly_the_hosts_it_names() -> None:
    institution = NetworkPolicy("institution", ("llm.lumc.nl", ".internal.lumc.nl"))
    for url in (
        "http://localhost:11434/v1",
        "https://llm.lumc.nl/v1",
        "https://a.internal.lumc.nl/x",
    ):
        institution.check(url)
    for url in (
        "https://api.openai.com/v1",
        "https://evil-lumc.nl/x",
        "https://internal.lumc.nl.evil.com/",
    ):
        with pytest.raises(ContractError, match="outside the institution"):
            institution.check(url)
    with pytest.raises(ContractError, match="offline"):
        NetworkPolicy().check("http://localhost:11434/v1")
    NetworkPolicy("online").check("https://api.openai.com/v1")
    with pytest.raises(ContractError, match="only http"):
        NetworkPolicy("online").check("file:///etc/passwd")


def test_an_allowed_host_cannot_redirect_a_request_outside() -> None:
    handler = _CheckedRedirect(NetworkPolicy("institution"))
    request = urllib.request.Request("http://localhost/x")
    with pytest.raises(ContractError, match="outside the institution"):
        handler.redirect_request(request, None, 302, "Found", {}, "https://example.com/y")


def test_offline_refuses_before_any_socket_is_opened(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("cell_curator.net._open", lambda *_: pytest.fail("opened a socket"))
    with pytest.raises(ContractError):
        urlopen(
            urllib.request.Request("https://example.com"), timeout=1, policy=NetworkPolicy()
        )


def test_offline_is_the_default_and_the_old_switch_still_loads() -> None:
    assert KnowledgeConfig().network == "offline"
    assert (
        KnowledgeConfig.model_validate({"allow_optional_remote_enhancement": True}).network
        == "online"
    )
    assert (
        KnowledgeConfig.model_validate({"allow_optional_remote_enhancement": False}).network
        == "offline"
    )
