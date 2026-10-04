"""
Builds an httpx.Client authenticated with Onshape's HMAC API-key scheme
(see auth.py), reading credentials from the repo root's .env (gitignored --
never commit it).

Uses plain httpx rather than onshape_rest_api_client.Client: the latter is
kept installed and useful as a reference for exact endpoint paths/request
shapes (it's generated straight from Onshape's own OpenAPI spec, under
onshape_rest_api_client/api/<tag>/<operation>.py), but its generated
response models have at least one real deserialization bug (crashes on a
null nested "company" field in BTUserOAuth2SummaryInfo, hit and confirmed
while validating auth here) -- not worth taking on that fragility for
every endpoint when a plain httpx.Client + manual response.json() works
and was already proven against the live API.
"""
import os
from pathlib import Path

import httpx
from dotenv import load_dotenv

from .auth import OnshapeApiKeyAuth

REPO_ROOT = Path(__file__).resolve().parent.parent


class _RateLimitedClient(httpx.Client):
    """Retries 429 Too Many Requests, honouring Retry-After when Onshape
    sends it and backing off exponentially when it does not.

    Added 2026-10-04 after a sweep lost 10 points in a row to 429s. Retrying
    here, at send(), rather than in a transport, means every attempt goes
    back through OnshapeApiKeyAuth.auth_flow and is re-signed with a fresh
    nonce and date -- replaying the first attempt's signed headers would be
    rejected as a reused nonce.
    """
    MAX_RETRIES = 8
    BASE_DELAY_S = 5.0
    MAX_DELAY_S = 300.0

    def send(self, request, **kwargs):
        import time
        for attempt in range(self.MAX_RETRIES + 1):
            resp = super().send(request, **kwargs)
            if resp.status_code != 429 or attempt == self.MAX_RETRIES:
                return resp
            try:
                delay = float(resp.headers.get("Retry-After", ""))
            except ValueError:
                delay = self.BASE_DELAY_S * 2 ** attempt
            if delay > self.MAX_DELAY_S:
                # A quota window, not a burst limit (seen 2026-10-04:
                # Retry-After 6975s with x-rate-limit-remaining 0). Retrying
                # inside it only burns time; fail now with the reset time.
                print(f"Onshape API quota exhausted on {request.method} {request.url.path}: "
                      f"Retry-After {delay:.0f}s (~{delay / 3600:.1f} h). Not retrying.",
                      flush=True)
                return resp
            delay = max(delay, 1.0)
            print(f"Onshape 429 on {request.method} {request.url.path}; retry "
                  f"{attempt + 1}/{self.MAX_RETRIES} in {delay:.0f}s", flush=True)
            resp.close()
            time.sleep(delay)


def get_client(timeout: float = 60.0) -> httpx.Client:
    load_dotenv(REPO_ROOT / ".env")
    return _RateLimitedClient(
        base_url=os.environ["ONSHAPE_BASE_URL"],
        auth=OnshapeApiKeyAuth(
            access_key=os.environ["ONSHAPE_ACCESS_KEY"],
            secret_key=os.environ["ONSHAPE_SECRET_KEY"],
        ),
        timeout=timeout,
    )
