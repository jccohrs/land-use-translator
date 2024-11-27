class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def print_section_heading(heading):
    print("______________________")
    print(heading)
    print("...")

def create_backgr_vars(nr_var, num_init):
    vars_list = ""
    for i in range(nr_var):
        num = i + num_init
        if len(str(num)) > 1:
           var_name = f"var8{num}"
        else:
            var_name = f"var80{num}"
        if i > 0:
            vars_list += "," + var_name
        else:
            vars_list += var_name
    return vars_list
