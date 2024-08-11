AUTH_PATH = '../Auth_file.txt'


class Auth():
    """
    Represents authentication information for a system.

    This class provides a simple way to store and retrieve authentication credentials,
    either from a specified text file or using default values.

    Attributes:
        username (str): The username for authentication.
        password (str): The password for authentication.

    Methods:
        from_txt(path: str) -> Auth:
            Creates an Auth instance by reading authentication information from a text file.

        from_default() -> Auth:
            Creates an Auth instance using default authentication information from the specified file path.

    Example:
        Creating an instance of Auth from a text file:

        >>> auth_instance = Auth.from_txt('../Auth_file.txt')

        Creating an instance of Auth using default values:

        >>> default_auth_instance = Auth.from_default()
    """

    def __init__(self, username, password):
        """
        Initializes an Auth instance with the provided username and password.

        Args:
            username (str): The username for authentication.
            password (str): The password for authentication.
        """
        self.username = username
        self.password = password

    def from_txt(path):
        """
        Creates an Auth instance by reading authentication information from a text file.

        Args:
            path (str): The path to the text file containing authentication information.

        Returns:
            Auth: An instance of the Auth class with credentials read from the text file.
        """
        f = open(path, "r")
        lines = f.readlines()
        for line in lines:
            if('username' in line):
                username = line.split('=')[1].replace('\n', '')
            if('password' in line):
                password = line.split('=')[1].replace('\n', '')
        return(Auth(username, password))

    def from_default():
        """
        Creates an Auth instance using default authentication information from the specified file path.

        Returns:
            Auth: An instance of the Auth class with credentials from the default file path.
        """
        return(Auth.from_txt(AUTH_PATH))
