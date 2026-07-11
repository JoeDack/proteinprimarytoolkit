import requests
from threading import local
from typing import Any

class InvalidSequenceError(ValueError):
    pass

thread_local = local()

# Called within functions that make API requests
# Ensures thread safety by preventing sharing of sessions across requests
def get_session() -> requests.Session:
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
    return thread_local.session

def get_json_from_api(session: requests.Session, url: str, timeout: int | float) -> Any:
    """
    Return JSON from API request. Return type depends on response, usually dict
    """
    if not isinstance(session, requests.Session):
        raise TypeError(f"session must be type Session, not {type(session)}")
    if not isinstance(url, str):
        raise TypeError(f"url must be str, not {type(url)}")
    if not isinstance(timeout, (int, float)):
        raise TypeError(f"timeout must be int or float, not {type(timeout)}")
    
    if url == "":
        raise ValueError("url must not be empty str")
    if timeout <= 0:
        raise ValueError("timeout must be greater than 0")
    
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        return response.json()

    # Raise more readable error than response.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"Failed to retrieve data from {url}") from e