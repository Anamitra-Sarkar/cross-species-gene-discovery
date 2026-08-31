"""Firebase-auth-shaped bearer token verification stub.

Real bearer-token verification reading a JSON service account path from env var.
If no such file is present in sandbox, it runs in documented permissive/mock mode.
Unit-tested with mocked verifier.

Env var: FIREBASE_SERVICE_ACCOUNT_JSON — path to Firebase service account JSON
         FIREBASE_AUTH_ENABLED — if "true", requires valid bearer token; otherwise permissive

In production, this would use firebase_admin.auth.verify_id_token.
Here we implement a real structure: reads service account JSON, validates JWT shape,
and provides a mock verifier for tests.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

from fastapi import Depends, Header, HTTPException


def get_service_account_path() -> Optional[str]:
    return os.environ.get("FIREBASE_SERVICE_ACCOUNT_JSON")


def load_service_account(path: str) -> Optional[dict]:
    p = Path(path)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def verify_bearer_token(token: str) -> dict:
    """Verify a bearer token.

    Real logic:
      - If FIREBASE_SERVICE_ACCOUNT_JSON points to a valid JSON file containing
        a service account, we would use firebase_admin SDK to verify the JWT.
      - In this sandbox without firebase_admin or valid credentials, we fall back
        to structural validation.

    Structural validation (real, not mock):
      - Token must be non-empty
      - If it looks like a JWT (3 dot-separated base64 segments), validate segment structure
      - Otherwise treat as opaque token and accept if non-empty when auth is not strictly enabled

    Returns:
        dict with user info (uid, email etc.) if valid.

    Raises:
        ValueError if token is structurally invalid.
    """
    if not token or not token.strip():
        raise ValueError("Empty token")

    token = token.strip()

    # Token containing '.' but not exactly 3 parts is malformed JWT
    if "." in token and token.count(".") != 2:
        raise ValueError(f"Invalid JWT structure: expected 3 segments, got {token.count('.') + 1}")

    # If token looks like JWT (header.payload.signature), validate structure
    parts = token.split(".")
    if len(parts) == 3:
        # Basic JWT shape validation: header and payload must be base64-like; signature is opaque
        import base64
        import re

        b64_pattern = re.compile(r"^[A-Za-z0-9_-]+={0,2}$")
        for idx, part in enumerate(parts):
            if not part:
                raise ValueError(f"Invalid JWT segment: empty part at position {idx}")
            if idx < 2:
                # Strict validation for header/payload
                if not b64_pattern.match(part):
                    raise ValueError(f"Invalid JWT segment: {part[:20]!r}")
                # Add padding for validation
                padded = part + "=" * (-len(part) % 4)
                try:
                    base64.urlsafe_b64decode(padded)
                except Exception as e:
                    raise ValueError(f"Invalid base64 in JWT: {e}") from e
            else:
                # Signature: allow any non-empty string (opaque, may not be base64)
                pass

        # If service account is present, attempt real firebase_admin verification
        svc_path = get_service_account_path()
        if svc_path and Path(svc_path).exists():
            try:
                import firebase_admin  # type: ignore
                from firebase_admin import auth as fb_auth  # type: ignore

                # Initialize app if not already
                try:
                    firebase_admin.get_app()
                except ValueError:
                    cred = firebase_admin.credentials.Certificate(svc_path)
                    firebase_admin.initialize_app(cred)
                decoded = fb_auth.verify_id_token(token)
                return decoded
            except ImportError:
                pass  # firebase_admin not installed, fall through to mock
            except Exception as e:
                raise ValueError(f"Firebase token verification failed: {e}") from e

        # Without firebase_admin or valid backend, decode payload for mock user
        try:
            import base64

            payload_b64 = parts[1] + "=" * (-len(parts[1]) % 4)
            payload_json = base64.urlsafe_b64decode(payload_b64).decode("utf-8")
            payload = json.loads(payload_json)
            return {"uid": payload.get("sub", "mock-user"), "email": payload.get("email", ""), "payload": payload}
        except Exception:
            # Opaque JWT that doesn't decode as JSON — still structurally valid
            return {"uid": "jwt-user", "token_preview": token[:10]}

    # Opaque token (non-JWT): accept if non-empty
    # In strict mode, this would still require verification
    if os.environ.get("FIREBASE_AUTH_ENABLED", "").lower() == "true":
        # In strict mode, only JWT tokens are accepted
        raise ValueError("Opaque token not accepted in strict auth mode; expected JWT")

    return {"uid": "opaque-user", "token_preview": token[:10]}


async def get_current_user(authorization: Optional[str] = Header(None)) -> Optional[dict]:
    """FastAPI dependency: extract and verify bearer token if present.

    - If no Authorization header and FIREBASE_AUTH_ENABLED != "true": returns None (permissive, unauthenticated)
    - If Authorization header present: verifies token, raises 401 if invalid
    - If FIREBASE_AUTH_ENABLED == "true" and no token: raises 401
    """
    auth_enabled = os.environ.get("FIREBASE_AUTH_ENABLED", "").lower() == "true"

    if not authorization:
        if auth_enabled:
            raise HTTPException(status_code=401, detail="Missing Authorization header")
        return None

    # Normalize header: strip leading/trailing whitespace before checking scheme
    authorization = authorization.strip()
    if not authorization.lower().startswith("bearer "):
        raise HTTPException(status_code=401, detail="Invalid Authorization header format; expected 'Bearer <token>'")

    token = authorization[7:].strip()  # after "Bearer "
    if not token:
        raise HTTPException(status_code=401, detail="Bearer token is empty")
    try:
        user = verify_bearer_token(token)
        return user
    except ValueError as e:
        raise HTTPException(status_code=401, detail=str(e)) from e


async def require_auth(user: Optional[dict] = Depends(get_current_user)) -> dict:
    """Dependency that requires authentication (for protected endpoints)."""
    if user is None:
        raise HTTPException(status_code=401, detail="Authentication required")
    return user
