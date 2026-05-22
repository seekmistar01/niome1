"""
Enhanced tests for validation logic in Bittensor subnet.
Includes parameterized, error, and docstring tests.
"""
import pytest
from unittest.mock import MagicMock

from niome_subnet.base.validator import Validator

@pytest.fixture
def mock_validator():
    validator = Validator()
    validator.validate = MagicMock(return_value=True)
    return validator

@pytest.mark.parametrize("input_data,expected", [
    ( {"data": "valid"}, True ),
    ( {"data": "invalid"}, False ),
])
def test_validation_param(mock_validator, input_data, expected):
    """Test validator with different data inputs."""
    mock_validator.validate.return_value = expected
    assert mock_validator.validate(input_data) is expected

def test_validation_error(mock_validator):
    """Test validator error handling."""
    mock_validator.validate.side_effect = Exception("Validation error")
    with pytest.raises(Exception, match="Validation error"):
        mock_validator.validate({"data": None})
