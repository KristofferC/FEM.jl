

function read_inp(filename)
    open(filename, "rU") do f

        # Read mesh objects
        num_elems = 0
        while True:
            l = readline(f)
            keyword = split(",")[1]

            if keyword == "*Part":
                mesh = _read_part(f)
            elif keyword == "*Node":
                _read_nodes(f, mesh)
            elif keyword == "*Element":
                num_elems += _read_elements(f, mesh, num_elems)
            elif keyword == "*Elset":
                _read_element_set(f, mesh)
            elif keyword == "*Nset":
                _read_node_set(f, mesh)
            elif keyword == "*End Part":
                break
            else:
                f.readline()
                continue

        return mesh


function _read_part(l):
    """
    Reads the part name and creates a mesh.
    """

    re_part = r"\*Part, name=(.*)"

    match = re_part.match(line)
    if not match:
        raise ReadInpFileError("Error parsing file. Expected '*Part, "
                               "name=XXX', read '" + line + "'.")

    return GeoMesh()


function _read_nodes(f, mesh):
    """
    Reads nodes from the file.
    """
    line = f.readline()
    if not (line == "*Node\n"):
        raise ReadInpFileError("\nError parsing file. Expected '*Node',"
                               " read '" + line + "'.")

    num_nodes = 0
    while True:
        start_of_line = f.tell()
        line = f.readline()
        if line.strip() == '':
            continue
        if line[0] == '*':
            f.seek(start_of_line)
            return
        num_nodes += 1
        if verbose == 1:
            print ("\rReading nodes, %d nodes read" % num_nodes),
        node_numbers = [to_number(x) for x in line.strip().split(',')]
        node = Node(node_numbers[0], node_numbers[1:])
        mesh.add_node(node)
        if verbose == 2:
            print ("Read {0}.\n".format(node))


function _read_elements(f, mesh, num_elems):
    """
    Reads elements from the file.
    """
    line = f.readline()
    re_element = re.compile("\*Element, type=(.*)")
    match = re_element.match(line)
    if not match:
        raise ReadInpFileError("\nError parsing file. Expected '*Element, \
        type=XXX', got '" + line + "'.")

    element_name = re_element.match(line).group(1)
    while True:
        start_of_line = f.tell()
        line = f.readline()
        if line.strip() == '':
            continue
        if line[0] == '*':
            f.seek(start_of_line)
            return num_elems
        num_elems += 1
        if verbose == 1:
            print ("\rReading element %s, with id %d."
                   % (element_name, num_elems)),

        element_numbers = map(to_number, line.strip().split(','))

        try:
            element_name_lolfem = element_translation[element_name]
        except KeyError:
            print ("ERROR: Unknown element: " + element_name + ", ignoring...")
        else:
            element_class = getattr(elements, element_name_lolfem)
            # Section number
            element = element_class(element_numbers[0], element_numbers[1:])
            mesh.add_element(element)


function _read_element_set(f, mesh):
    """
    Reads element sets from the file.
    """
    line = f.readline()
    re_element_set = re.compile("\*Elset, elset=(.*)")
    match = re_element_set.match(line)
    if not match:
        raise ReadInpFileError("Error parsing file. Expected '*Elset, "
                               "elset=X', got '" + line + "'.")

    element_set_name = re_element_set.match(line).group(1)

    if element_set_name.startswith("face"):
        dim = 2
    elif element_set_name.startswith("poly"):
        dim = 3
    else:
        dim = None
    if verbose == 1 or verbose == 2:
        print ("\rReading element set {0:s}.".format(element_set_name)),

    full_str = ""
    if element_set_name.endswith("generate"):
        element_set_name = element_set_name[0:-10]
        element_set = ElementSet(element_set_name, dim)
        line = f.readline().strip()
        generate_info = map(to_number, line.split(','))
        start, stop, step = generate_info[
            0], generate_info[1], generate_info[2]
        element_set.ids = range(start, stop + 1, step)
        mesh.element_sets[element_set_name] = element_set
        return
    else:
        element_set = ElementSet(element_set_name, dim)
        while True:
            start_of_line = f.tell()
            line = f.readline()
            if line.strip() == '':
                continue
            if line[0] == '*':
                element_list = full_str.split(',')
                element_list = [item for item in element_list if item]
                element_set.ids = map(to_number, element_list)
                mesh.element_sets[element_set_name] = element_set
                f.seek(start_of_line)
                return
                # Read element ids until empty line
            full_str += line.strip() + ","


function _read_node_set(f, mesh):
    """
    Reads node sets from the file.
    """
    line = f.readline()
    re_node_set = re.compile("\*Nset, nset=(.*)")
    match = re_node_set.match(line)
    if not match:
        raise ReadInpFileError("Error parsing file. Expected '*Nset, "
                               "nset=X', got '" + line + "'.")
    node_set_name = re_node_set.match(line).group(1)
    if verbose == 1 or verbose == 2:
        print ("\rReading node set {0:s}.".format(node_set_name)),
    full_str = ""
    if node_set_name.endswith("generate"):
        node_set_name = node_set_name[0:-10]
        node_set = NodeSet(node_set_name)
        line = f.readline().strip()
        generate_info = map(to_number, line.split(','))
        start, stop, step = generate_info[
            0], generate_info[1], generate_info[2]
        node_set.ids = range(start, stop + 1, step)
        mesh.node_sets[node_set_name] = node_set
        return
    else:
        node_set = NodeSet(node_set_name)
        while True:
            start_of_line = f.tell()
            line = f.readline()
            if line.strip() == '':
                continue
            if line[0] == '*':
                # Remove empty strings
                node_list = full_str.split(',')
                # Remove empty strings
                node_list = [item for item in node_list if item]
                node_set.ids = map(to_number, node_list)
                mesh.node_sets[node_set_name] = node_set
                f.seek(start_of_line)
                return
            full_str += line.strip() + ","
