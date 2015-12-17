from validate_email import validate_email

from werkzeug.utils import secure_filename
from werkzeug.datastructures import FileStorage

from server.webapp.exceptions import InvalidFileType
from server.webapp.utils import allowed_file

def email(email_str):
    if validate_email(email_str):
        return email_str

    raise ValueError('{} is not a valid email'.format(email_str))


def hashed_password(password_str):
    if isinstance(password_str, basestring) and len(password_str) == 56:
        return password_str

    raise ValueError('Invalid password - expecting SHA224 - Received {} of length {} and type {}'.format(password_str, len(password_str), type(password_str)))


def secure_filename_input(orig_name):
    return secure_filename(orig_name)


def Json(orig):
    return orig


class SubRequest:

    def __init__(self, orig_dict):
        self.json = orig_dict


class SubParser:

    __name__ = 'Nested parser'

    def __init__(self, child_parser):
        self.child_parser = child_parser

    def __call__(self, item_to_parse):
        return self.child_parser.parse_args(req=SubRequest(item_to_parse))


class AllowedFileTypeMixin(object):
    "Mixin used of FileStorage subclasses to check the uploaded filetype"

    def __init__(self, *args, **kwargs):
        super(AllowedFileTypeMixin, self).__init__(*args, **kwargs)
        if not allowed_file(self.filename):
            raise ValueError('File type of {} is not accepted!'.format(self.filename))


class SafeFilenameStorage(FileStorage):

    def __init__(self, *args, **kwargs):
        super(SafeFilenameStorage, self).__init__(*args, **kwargs)
        if self.filename == 'file' and hasattr(self.stream, 'filename'):
            self.filename = self.stream.filename
        self.source_filename = self.filename
        self.filename = secure_filename(self.filename)


class AllowedFiletypeStorage(AllowedFileTypeMixin, FileStorage):
    pass


class AllowedSafeFilenameStorage(AllowedFileTypeMixin, SafeFilenameStorage):
    pass
