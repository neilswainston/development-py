'''
Unique (c) University of Manchester 2015

Unique is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import uuid

from flask import Flask

# Configuration:
DEBUG = True
SECRET_KEY = str(uuid.uuid4())

# Create application:
_APP = Flask(__name__)
_APP.config.from_object(__name__)


@_APP.route('/unique/<number>')
def unique(number):
    '''Unique identifier generator.'''
    return '\n'.join([str(uuid.uuid4()) for _ in range(0, int(number))])


if __name__ == '__main__':
    _APP.run(threaded=True)
