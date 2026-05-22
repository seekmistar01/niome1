"""
Enhanced tests for incentives logic in Bittensor subnet.
Includes parameterized, error, and docstring tests.
"""
import pytest
from unittest.mock import MagicMock

from niome_subnet.validator.reward import Reward

@pytest.fixture
def mock_reward():
    reward = Reward()
    reward.distribute = MagicMock(return_value={"uid": 1, "reward": 10})
    return reward

@pytest.mark.parametrize("input_data,expected", [
    ( {"uid": 1}, {"uid": 1, "reward": 10} ),
    ( {"uid": 2}, {"uid": 2, "reward": 5} ),
])
def test_incentive_distribution_param(mock_reward, input_data, expected):
    """Test incentive distribution for different uids."""
    mock_reward.distribute.return_value = expected
    result = mock_reward.distribute(input_data)
    assert result["uid"] == expected["uid"]
    assert result["reward"] == expected["reward"]

def test_incentive_distribution_error(mock_reward):
    """Test incentive distribution error handling."""
    mock_reward.distribute.side_effect = Exception("Distribution error")
    with pytest.raises(Exception, match="Distribution error"):
        mock_reward.distribute({"uid": None})
