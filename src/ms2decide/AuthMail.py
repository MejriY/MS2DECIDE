AUTH_PATH = '../Auth_file.txt'


class AuthMail():
    """
    Represents authentication information for a mail system.

    This class provides a way to store and retrieve authentication credentials,
    including username, password, and mail address, either from a specified text file or using default values.

    Attributes:
        username (str): The username for authentication.
        password (str): The password for authentication.
        mail (str): The email address for authentication.

    Methods:
        from_txt(path: str) -> AuthMail:
            Creates an AuthMail instance by reading authentication information from a text file.

        from_default() -> AuthMail:
            Creates an AuthMail instance using default authentication information from the specified file path.

    Example:
        Creating an instance of AuthMail from a text file:

        >>> auth_mail_instance = AuthMail.from_txt('../Auth_file.txt')

        Creating an instance of AuthMail using default values:

        >>> default_auth_mail_instance = AuthMail.from_default()
    """

    def __init__(self, username, password, mail):
        """
        Initializes an AuthMail instance with the provided username, password, and mail.

        Args:
            username (str): The username for authentication.
            password (str): The password for authentication.
            mail (str): The email address for authentication.
        """
        self.username = username
        self.password = password
        self.mail = mail

    def from_txt(path):
        """
        Creates an AuthMail instance by reading authentication information from a text file.

        Args:
            path (str): The path to the text file containing authentication information.

        Returns:
            AuthMail: An instance of the AuthMail class with credentials read from the text file.
        """
        f = open(path, "r")
        lines = f.readlines()
        for line in lines:
            if('username' in line):
                username = line.split('=')[1].replace('\n', '')
            if('password' in line):
                password = line.split('=')[1].replace('\n', '')
            if('mail' in line):
                mail = line.split('=')[1].replace('\n', '')
        return(AuthMail(username, password, mail))

    def from_default():
        """
        Creates an AuthMail instance using default authentication information from the specified file path.

        Returns:
            AuthMail: An instance of the AuthMail class with credentials from the default file path.
        """
        return(AuthMail.from_txt(AUTH_PATH))
