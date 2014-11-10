def legitimize(filename='vasprun.xml'):
    """Add closing tags </calculation> and </modeling> to the tree when job is walled.
    Takes a filename.
    """
    with open(filename, 'r+') as f:
        lines = f.readlines()
        if lines[-1] == '</scstep>':
            lines.pop(-1)
        lines.append('</calculation>\n')
        lines.append('</modeling>\n')
        f.seek(0)
        f.write("".join(lines))


def parse(filename='vasprun.xml'):
    """
    Takes a filename.
    Returns the root element.
    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(filename)
    return tree.getroot()


def iterprint(elem, xpath='.'):
    """
    Takes an xml element, an xpath string.
    Prints element's info in a more straightforward way.
    """
    for e in elem.findall(xpath):
        print e.tag, e.attrib,
        if e.text:
            text = e.text.strip()
            print text
        else:
            print