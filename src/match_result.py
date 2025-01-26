from dataclasses import dataclass

@dataclass
class MatchResult:
    score: float = 0.0
    read: str = ''
    # primer: str
    start_index: int = -1
    end_index: int = -1

