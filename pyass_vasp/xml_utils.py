def legitimize_vasprun_xml(filename):
    """Add closing tags </calculation> and </modeling> to the tree when job is walled.

    Takes a filename.
    """
    with open(filename, 'r+') as f:
        lines = f.readlines()
        lines.pop(-1)
        lines.append('</calculation>\n')
        lines.append('</modeling>\n')
        f.seek(0)
        f.write("".join(lines))

def parse(filename):
    """
    Takes a filename, returns the root.
    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(filename)
    return tree.getroot()

def iterprint(elem, xpath='.'):
    """
    Takes an xml element, an xpath string.
    Prints element's text.
    """
    for e in elem.findall(xpath):
        print "{0}: {1}".format(e.tag, e.text)