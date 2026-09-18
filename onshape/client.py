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


def get_client(timeout: float = 60.0) -> httpx.Client:
    load_dotenv(REPO_ROOT / ".env")
    return httpx.Client(
        base_url=os.environ["ONSHAPE_BASE_URL"],
        auth=OnshapeApiKeyAuth(
            access_key=os.environ["ONSHAPE_ACCESS_KEY"],
            secret_key=os.environ["ONSHAPE_SECRET_KEY"],
        ),
        timeout=timeout,
    )
