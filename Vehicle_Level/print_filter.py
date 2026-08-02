import builtins
import contextlib

YELLOW = '\033[93m'
RED = '\033[91m'
BOLD = '\033[1m'
ENDC = '\033[0m'
GREEN = '\033[92m'

def Color_MoS_Text(MoS_value, full_string):
    if MoS_value >= 0:
        text_color = GREEN
    else:
        text_color = RED
        
    colored_full_string = f"{text_color}{full_string}{ENDC}"
        
    return(colored_full_string)

@contextlib.contextmanager
def context_manager(print_everything: bool = True, print_margins: bool = True, print_titles: bool = True):
    original_print_function = builtins.print

    def filtered_print(*arguments, **keyword_arguments):
        is_margin = bool(keyword_arguments.pop("i_am_a_margin", False))
        is_title = bool(keyword_arguments.pop("i_am_a_title", False))
        always_print_this = bool(keyword_arguments.pop("always_print_this", False))
        
        if print_everything or always_print_this:
            return original_print_function(*arguments, **keyword_arguments)
        else:
            if (print_margins and is_margin):
                return original_print_function(*arguments, **keyword_arguments)
            elif (print_titles and is_title):
                return original_print_function(*arguments, **keyword_arguments)

    builtins.print = filtered_print
    try:
        yield
    finally:
        builtins.print = original_print_function
