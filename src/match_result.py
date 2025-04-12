from dataclasses import dataclass

@dataclass
class MatchResult:
    score: float = 0.0
    read: str = ''
    start_index: int = -1
    end_index: int = -1

    def is_mismatch(self):
        return True if self.start_index == -1 else False

