
def main():
    source(findFile("scripts", "test_helpers.py"))
    source(findFile("scripts", "common_checks.py"))
    startApplication("VatesSimpleGui")
    openFile("")
    check_mode_buttons(std=False, ms=True, ts=True, sp=True)
    check_time_controls(True)
    quitProgram()

