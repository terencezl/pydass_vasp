def legitimize_vasprun_xml(filename):
    """
    Takes a filename, makes vasprun.xml legitimate.
    """
    with open(filename, 'r+') as f:
        lines = f.readlines()
        lines.pop(-1)
        lines.append('</calculation>\n')
        lines.append('</modeling>\n')
        f.seek(0)
        f.write("".join(lines))


def parse_xml(filename):
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