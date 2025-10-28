class ProgressTracker:
    """Tracks and visualizes construction progress."""

    def __init__(self):
        self.progress_data = {}
        self.current_step = 0
        self.total_steps = 0

    def start_construction(self, total_steps):
        """Start tracking construction with total steps."""
        self.total_steps = total_steps
        self.current_step = 0
        self.progress_data.clear()

    def update_progress(self, step_name, step_data=None):
        """Update construction progress."""
        self.current_step += 1
        self.progress_data[step_name] = {
            'step_number': self.current_step,
            'data': step_data,
            'progress': (self.current_step / self.total_steps) * 100
        }

    def get_progress(self):
        """Get current progress percentage."""
        if self.total_steps == 0:
            return 0
        return (self.current_step / self.total_steps) * 100

    def generate_progress_report(self):
        """Generate progress report."""
        return {
            'current_step': self.current_step,
            'total_steps': self.total_steps,
            'progress_percentage': self.get_progress(),
            'step_details': self.progress_data
        }