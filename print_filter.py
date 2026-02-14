import builtins
import contextlib

@contextlib.contextmanager
def context_manager(print_everything: bool = True, print_margins: bool = True, print_titles: bool = True):
    original_print_function = builtins.print

    def filtered_print(*arguments, **keyword_arguments):
        is_margin = bool(keyword_arguments.pop("i_am_a_margin", False))
        is_title = bool(keyword_arguments.pop("i_am_a_title", False))
        
        if print_everything:
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
