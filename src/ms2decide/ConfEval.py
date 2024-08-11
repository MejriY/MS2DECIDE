from .Util import owa, Weights
from decimal import Decimal


class ConfEval():
    """
    Estimator class for calculating confidence and recommendation scores.

    Args:
        similarities (dict): Dictionary of similarity values for each tool.
        tanimotos (Tanimotos): Instance of Tanimotos class.
        a_priori_confidences (dict): A priori confidence values for each tool.
        weight_tanimoto (float): Weight for Tanimoto in the overall confidence computation.
        weights (list): List of weights for each tool.
        default_cos (float): Default cosine value.

    Attributes:
        similarities (dict): Dictionary of similarity values for each tool.
        tanimotos (Tanimotos): Instance of Tanimotos class.
        a_priori_confidences (dict): A priori confidence values for each tool.
        weight_tanimoto (float): Weight for Tanimoto in the overall confidence computation.
        weights (list): List of weights for each tool.
        default_cos (float): Default cosine value.
        tools (list): List of tools.
        relative_confidences (dict): Dictionary of relative confidences for each tool.
        confidences (dict): Dictionary of confidences for each tool.

    Raises:
        Exception: If similarity dict and tanimoto dict do not have the same keys.

    Methods:
        _confidences: Calculate confidences for each tool.
        _relative_confidences: Calculate relative confidences for each tool.
        confidence: Get the confidence score for a specific tool.
        relative_confidence: Get the relative confidence score for a specific tool.
        f: Calculate the f score.
        g: Calculate the g score.
        recommendation: Get the recommendation score.

    Example:
        # Create an Estimator instance
        estimator = Estimator(similarities, tanimotos, a_priori_confidences, weight_tanimoto, weights, default_cos)

        # Get confidence for a specific tool
        confidence_score = estimator.confidence("Tool1")

        # Get relative confidence for a specific tool
        relative_confidence_score = estimator.relative_confidence("Tool2")

        # Get the recommendation score
        recommendation_score = estimator.recommendation()
    """

    def __init__(self, similarities, tanimotos):
        """
        Initialize the Estimator.

        Args:
            similarities (dict): Dictionary of similarity values for each tool.
            tanimotos (Tanimotos): Instance of Tanimotos class.
            a_priori_confidences (dict): A priori confidence values for each tool.
            weight_tanimoto (float): Weight for Tanimoto in the overall confidence computation.
            weights (list): List of weights for each tool.
            default_cos (float): Default cosine value.
        """
        a_priori_confidences, weight_tanimoto, weights, default_cos = Weights()
        if(set(similarities.keys()) == set(tanimotos.data.keys())):
            self.similarities = similarities
            self.tanimotos = tanimotos
            self.a_priori_confidences = a_priori_confidences
            self.weight_tanimoto = weight_tanimoto
            self.weights = weights
            self.default_cos = default_cos
            self.tools = list(similarities.keys())
            self.relative_confidences = None
            self.confidences = None
        else:
            raise Exception(
                'similarity dict and tanimoto dict should have the same keys')

    def _confidences(self):
        """
        Calculate and set the confidences for each tool.
        """
        tools = self.tools
        similarities = self.similarities
        for i in similarities:
            if(str(similarities[i]) == 'nan'):
                similarities[i] = 0
        w_t = self.weight_tanimoto
        d = {i: 0 for i in tools}
        for tool in d:
            a = self.tanimotos.average_tool(tool)
            d[tool] = float(Decimal((1-w_t)*similarities[tool]+w_t * a))
        self.confidences = d

    def _relative_confidences(self):
        """
        Calculate and set the relative confidences for each tool.
        """
        if(self.confidences == None):
            self._confidences()
        confidences = self.confidences
        st = self.a_priori_confidences
        tools = self.tools
        s = float(Decimal(sum([Decimal(confidences[i]*st[i]) for i in tools])))
        d = {}
        for i in tools:
            if(s!=0):
                d[i]= float(Decimal((confidences[i] * st[i])/s))
            else:
                d[i]=0
        self.relative_confidences = d

    def confidence(self, tool):
        """
        Get the confidence value for a specific tool.

        Args:
            tool (str): Tool name.

        Returns:
            float: Confidence value for the tool.
        """
        if(self.confidences == None):
            self._confidences()
        return(self.confidences[tool])

    def relative_confidence(self, tool):
        """
        Get the relative confidence value for a specific tool.

        Args:
            tool (str): Tool name.

        Returns:
            float: Relative confidence value for the tool.
        """
        if(self.relative_confidences == None):
            self._relative_confidences()
        return(self.relative_confidences[tool])

    def f(self):
        """
        Compute and return the existence of our given spectrum in all tools.

        Returns:
            float: The existence of our given spectrum in all tools score.
        """
        if(len(self.tools) == 1):
            c = self.similarities[self.tools[0]]
            return(float(Decimal(c*c+(1-c)*self.default_cos)))
        elif(set(self.similarities.values()) == {0}) and (set(self.tanimotos.data.values()) == {0}):
            return(self.default_cos)
        else:
            tools = self.tools
            if(self.confidences == None):
                self._confidences()
            similarities = self.similarities
            if(self.relative_confidences == None):
                self._relative_confidences()
            rt = self.relative_confidences
            c_m = self.default_cos
            c = float(
                Decimal(sum(Decimal(similarities[tool]*rt[tool]) for tool in tools)))
            g = self.g()
# TODO il n’est pas approprié de renvoyer des résultats arrondis, c’est à l’appelant de choisir. Concernant le problème qui est j’imagine à l’origine de ce changement, lire par exemple https://davidamos.dev/the-right-way-to-compare-floats-in-python/
            return(round(g*c+(1-g)*c_m, 6))

    def g(self):
        """
        Compute and return the confidence of tools prediction.

        Returns:
            float: The confidence of tools prediction.
        """
        if(len(self.tools) == 1):
            return(self.similarities[self.tools[0]])
        elif(set(self.similarities.values()) == {0}) and (set(self.tanimotos.data.values()) == {0}):
            return(0)
        else:
            if(self.confidences == None):
                self._confidences()
            zetas = self.confidences
            w_g = self.weights
            return(round(owa(w_g, zetas.values()), 6))

    def recommendation(self):
        """
        Compute and return the recommendation score.

        Returns:
            float: The recommendation score.
        """
        g = self.g()
        f = self.f()
        return(round((1-f)+(1-g)*0.3, 6))
