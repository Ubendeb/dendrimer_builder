"""
Main entry point for Dendrimer Builder application.
"""

from src import DendrimerBuilder, FragmentInterface, ConstructionInterface


def main():
    """Main application function."""
    print("Dendrimer Builder System")
    print("=" * 30)

    # Initialize components
    builder = DendrimerBuilder()
    fragment_interface = FragmentInterface()
    construction_interface = ConstructionInterface(builder)

    print("System initialized successfully!")

    # Example usage pattern
    print("\nAvailable modules:")
    print("- Fragment management")
    print("- Dendrimer construction")
    print("- Structure visualization")
    print("- Data import/export")


if __name__ == "__main__":
    main()