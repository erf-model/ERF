"""Custom Sphinx directive to include remote files via URL"""
from docutils import nodes
from docutils.parsers.rst import Directive, directives
import requests

class RemoteInclude(Directive):
    required_arguments = 1  # URL
    optional_arguments = 0
    option_spec = {
        'language': directives.unchanged,
    }

    def run(self):
        url = self.arguments[0]
        language = self.options.get('language', 'text')

        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            content = response.text
        except Exception as e:
            error_msg = f"Failed to fetch {url}: {str(e)}"
            return [nodes.error(None, nodes.paragraph(text=error_msg))]

        literal = nodes.literal_block(content, content)
        literal['language'] = language
        return [literal]

def setup(app):
    app.add_directive('remote-include', RemoteInclude)
    return {'version': '0.1', 'parallel_read_safe': True}
