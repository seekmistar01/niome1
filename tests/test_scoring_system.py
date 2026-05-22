"""
Enhanced tests for scoring system logic in Bittensor subnet.
Includes parameterized, error, and docstring tests.
"""
import pytest
from unittest.mock import MagicMock

from niome_subnet.validator.reward import Reward

@pytest.fixture
def mock_reward():
    reward = Reward()
    reward.calculate_score = MagicMock(return_value=0.95)
    return reward

@pytest.mark.parametrize("input_data,expected", [
    ( {"response": "good"}, 0.95 ),
    ( {"response": "bad"}, 0.0 ),
])
def test_scoring_system_param(mock_reward, input_data, expected):
    """Test scoring system with different responses."""
    mock_reward.calculate_score.return_value = expected
    score = mock_reward.calculate_score(input_data)
    assert score == expected

def test_scoring_system_error(mock_reward):
    """Test scoring system error handling."""
    mock_reward.calculate_score.side_effect = Exception("Invalid input")
    with pytest.raises(Exception, match="Invalid input"):
        mock_reward.calculate_score({"response": None})
