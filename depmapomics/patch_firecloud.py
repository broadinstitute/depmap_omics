import random
import time
import requests
import firecloud.api
from functools import update_wrapper

def wrap_call(obj, prop, wrapper):    
    fn = getattr(obj, prop)
    call = lambda *args, **kwargs: wrapper(fn, args, kwargs)
    update_wrapper(wrapper, call)
    setattr(obj, prop, call)

def random_connection_fault_wrapper(probability_of_failure=0.5):
    """
    A wrapper which can be installed via install_patches. This one will randomly throw 
    ConnectionErrors to simulate failures.
    """
    def call_fn_with_random_errors(fn, args, kwargs):
        if random.random() < probability_of_failure:
            raise requests.ConnectionError()
        return fn(*args, **kwargs)
    
    return call_fn_with_random_errors

    
def retry_on_connection_error_wrapper(max_retries=10, initial_sleep=0.1):
    """
    Creates a wrapper which can be installed via install_patches. This one will retry the call 
    (with exponential backoff) if an ConnectionError is thrown
    """
    def call_fn_with_retries(fn, args, kwargs):
        sleep_time = initial_sleep
        attempt = 0
        while True:
            try:
                return fn(*args, **kwargs)
            except requests.ConnectionError as ex:
                attempt += 1
                if attempt >= max_retries:
                    break
                print(f"Got exception {ex}, retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
                sleep_time *= 2
        raise Exception("Too many ConnectionErrors, giving up.")
    return call_fn_with_retries

def install_patches(patches=[retry_on_connection_error_wrapper()]):
    """
    Installs one or more wrappers on the "get" method of the AuthorizedSession instance that 
    FISS uses to make firecloud api calls.
    """
    assert firecloud.api.__SESSION is None, "firecloud.api.__SESSION has already been set. Calls to install_patches() must be called before any firecloud api calls."
    
    def monkeypatch_get_method(authorized_session_fn, args, kwargs):
        session = authorized_session_fn(*args, **kwargs)
        for patch in patches:
            wrap_call(session, "get", patch)
        return session
    
    wrap_call(firecloud.api, "AuthorizedSession", monkeypatch_get_method)

