# tests/test_attributes.py
import pytest
from dephycf import attributes

def test_required_in_known():
    """
    Verify that all required attributes are present in the list of known attributes.
    """
    missing = [attr for attr in attributes.required_attributes 
               if attr not in attributes.known_attributes]
    assert not missing, f"The following required attributes are missing from known_attributes: {missing}"

def test_no_duplicates_known():
    """
    Verify that there are no duplicate entries in known_attributes.
    """
    duplicates = set([x for x in attributes.known_attributes if attributes.known_attributes.count(x) > 1])
    assert not duplicates, f"Duplicate entries found in known_attributes: {duplicates}"

def test_all_strings():
    """
    Verify that all entries in known_attributes and required_attributes are strings.
    """
    non_strings_known = [x for x in attributes.known_attributes if not isinstance(x, str)]
    non_strings_required = [x for x in attributes.required_attributes if not isinstance(x, str)]

    assert not non_strings_known, f"Non-string entries in known_attributes: {non_strings_known}"
    assert not non_strings_required, f"Non-string entries in required_attributes: {non_strings_required}"

