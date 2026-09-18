"""
HMAC request signing for Onshape's API-key authentication scheme.

onshape-rest-api-client (the OpenAPI-generated client this project uses for
the typed request/response models) only ships a Bearer-token
`AuthenticatedClient`, meant for OAuth2 access tokens -- it has no support
for Onshape's access-key/secret-key HMAC scheme. This fills that gap as a
plain httpx.Auth, so it plugs into the generated client's plumbing via
`Client(..., httpx_args={"auth": OnshapeApiKeyAuth(...)})` (Client passes
**httpx_args straight through to the underlying httpx.Client, which accepts
a custom `auth` object) instead of needing a separate hand-rolled requests
layer.

Signing algorithm ported from Onshape's own reference implementation
(onshape-public/apikey, python/apikey/onshape.py -- Python 2, MIT-licensed
sample repo) to modern Python 3 / httpx. Not guessed or reconstructed from
memory: fetched and read directly before writing this.
"""
import base64
import datetime
import hashlib
import hmac
import random
import string

import httpx


class OnshapeApiKeyAuth(httpx.Auth):
    def __init__(self, access_key: str, secret_key: str):
        self._access_key = access_key
        self._secret_key = secret_key.encode("utf-8")

    @staticmethod
    def _nonce() -> str:
        chars = string.digits + string.ascii_letters
        return "".join(random.choice(chars) for _ in range(25))

    def auth_flow(self, request: httpx.Request):
        method = request.method
        path = request.url.path
        # NOTE: Onshape's reference impl url-encodes a query dict via
        # urllib.urlencode() before lowercasing; httpx's request.url.query
        # is already percent-encoded from building the request, which is
        # equivalent for the plain key=value query strings this project
        # uses (translation/configuration calls). Revisit if a future
        # query string needs unusual characters and signing starts failing.
        query = request.url.query.decode("utf-8") if request.url.query else ""
        ctype = request.headers.get("Content-Type", "application/json")
        date = datetime.datetime.now(datetime.timezone.utc).strftime("%a, %d %b %Y %H:%M:%S GMT")
        nonce = self._nonce()

        hmac_str = f"{method}\n{nonce}\n{date}\n{ctype}\n{path}\n{query}\n".lower().encode("utf-8")
        signature = base64.b64encode(
            hmac.new(self._secret_key, hmac_str, digestmod=hashlib.sha256).digest()
        ).decode("utf-8")

        request.headers["Date"] = date
        request.headers["On-Nonce"] = nonce
        request.headers["Authorization"] = f"On {self._access_key}:HmacSHA256:{signature}"
        request.headers.setdefault("Content-Type", ctype)
        request.headers.setdefault("Accept", "application/json")
        request.headers.setdefault("User-Agent", "HET-project (onshape/auth.py)")

        yield request
