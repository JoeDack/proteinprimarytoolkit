import requests
from threading import local
from typing import Any

DEFAULT_TIMEOUT = 30

class InvalidSequenceError(ValueError):
    pass

class APIError(RuntimeError):
    pass

class InvalidPDBIDError(ValueError):
    pass

_thread_local = local()

def raise_for_api_status(response: requests.Response) -> None:
    try:
        response.raise_for_status()
    # Raise more readable error than response.raise_for_status()
    except requests.HTTPError as e:
        raise APIError(f"HTTP {response.status_code} ({response.reason}) returned from {response.url}") from e


# Called within functions that make API requests
# Ensures thread safety by preventing sharing of sessions across requests
def get_session() -> requests.Session:
    if not hasattr(_thread_local, "session"):
        _thread_local.session = requests.Session()
    return _thread_local.session

def get_response(session: requests.Session, url: str, timeout: int | float = DEFAULT_TIMEOUT) -> requests.Response:
    if not isinstance(url, str):
        raise TypeError(f"url must be str, not {type(url)}")
    # Checking if bool because bool is a subclass of int in Python
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise TypeError(f"timeout must be int or float, not {type(timeout)}")
    
    if not url.strip():
        raise ValueError("url must not be blank")
    if timeout <= 0:
        raise ValueError("timeout must be greater than 0")
    
    try:
        response = session.get(url, timeout=timeout)
        raise_for_api_status(response)
        return response

    except requests.Timeout as e:
        raise APIError(f"Request timed out: {url}") from e
    except requests.ConnectionError as e:
        raise APIError(f"Failed to connect: {url}") from e
    except requests.RequestException as e:
        raise APIError(f"Request failed: {url}") from e
    
def get_json_from_api(session: requests.Session, url: str, timeout: int | float = DEFAULT_TIMEOUT) -> Any:
    """
    Return JSON from API request. Return type depends on response, usually dict
    """
    response = get_response(session, url, timeout)
    try:
        return response.json()
    except requests.exceptions.JSONDecodeError as e:
        raise APIError(f"Invalid JSON returned from {url}") from e

def get_text_from_api(session: requests.Session, url: str, timeout: int | float = DEFAULT_TIMEOUT) -> str:
    response = get_response(session, url, timeout)
    return response.text