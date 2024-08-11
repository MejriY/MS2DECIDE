class MatchedSpectra:
    """
    Represents a matched spectrum along with associated information.

    This class is used to store information about a matched spectrum, including the spectrum itself,
    the InChI (International Chemical Identifier) of the corresponding compound, and the matching score.

    Attributes:
        spectrum: The matched spectrum data.
        inchi (str): The InChI (International Chemical Identifier) of the corresponding compound.
        score (float): The matching score indicating the degree of similarity between the spectrum and the compound.

    Example:
        Creating an instance of MatchedSpectra:

        >>> spectrum_data = [...]  # Replace with actual spectrum data
        >>> inchi_code = "InChI=1S/C6H12O6/c7-1-3-4-5(9)6(10)8-2-1/h1,3-6,9-10H,2H2"  # Replace with actual InChI code
        >>> matching_score = 0.85  # Replace with actual matching score
        >>> matched_spectrum = MatchedSpectra(spectrum_data, inchi_code, matching_score)
    """

    def __init__(self, sp, inchi, score):
        """
        Initializes a MatchedSpectra instance.

        Args:
            spectrum: The matched spectrum data.
            inchi (str): The InChI (International Chemical Identifier) of the corresponding compound.
            score (float): The matching score indicating the degree of similarity between the spectrum and the compound.
        """
        self.spectrum = sp
        self.inchi = inchi
        self.score = score
