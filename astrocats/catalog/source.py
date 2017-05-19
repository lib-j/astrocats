"""Class for representing sources of data.
"""
import re

import arxiv

from astrocats.catalog import utils
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class SOURCE(KeyCollection):
    """`KeyCollection` for the `Source` class.

    Attributes
    ----------
    NAME : STRING
    BIBCODE : STRING
    URL : STRING
    ACKNOWLEDGMENT : STRING
    REFERENCE : STRING
    ALIAS : NUMERIC
        Numerical alias (shorthand) for this entry.  Saved as a string (or
        list of strings), despite being stored as an integer.
    SECONDARY : BOOL
        Whether the given source is one which collected data from another,
        'Primary'-source, from which it actually originated

    """
    # Strings
    NAME = Key('name', KEY_TYPES.STRING)
    BIBCODE = Key('bibcode', KEY_TYPES.STRING)
    ARXIVID = Key('arxivid', KEY_TYPES.STRING)
    DOI = Key('doi', KEY_TYPES.STRING)
    URL = Key('url', KEY_TYPES.STRING, compare=False)
    ACKNOWLEDGMENT = Key('acknowledgment', KEY_TYPES.STRING, compare=False)
    REFERENCE = Key('reference', KEY_TYPES.STRING, compare=False)
    # Numbers
    ALIAS = Key('alias', KEY_TYPES.INTEGER, compare=False)
    # Booleans
    SECONDARY = Key('secondary', KEY_TYPES.BOOL, compare=False)
    PRIVATE = Key('private', KEY_TYPES.BOOL, compare=False)


class Source(CatDict):
    """Representation for the source/attribution of a data element.
    """

    _KEYS = SOURCE
    _MAX_AUTHORS_TO_LIST = 3

    def __init__(self, parent, derive_parameters=False, **kwargs):
        self._REQ_KEY_SETS = [
            [SOURCE.ALIAS],
            [SOURCE.BIBCODE, SOURCE.ARXIVID, SOURCE.URL, SOURCE.NAME]
        ]
        super(Source, self).__init__(parent, **kwargs)
        if derive_parameters:
            self.ads = utils.import_ads()
            self.derive_parameters()

        return

    def derive_parameters(self):
        self._log.debug("Source.derive_parameters()")
        bib = self.get(SOURCE.BIBCODE)
        arx = self.get(SOURCE.ARXIVID)
        doi = self.get(SOURCE.DOI)
        url = self.get(SOURCE.URL)
        name = self.get(SOURCE.NAME)
        authors = None
        year = None
        self._log.info("name = '{}', bib = '{}', arx = '{}', url = '{}', doi = '{}'".format(
            name, bib, arx, url, doi))

        if None not in [bib, arx, doi, url, name]:
            self._log.debug("Source entry is already complete.")
            return

        # Use the URL to try to get a bibcode or arxiv ID
        #    determine if URL is useful here
        from_url = (bib is None) or (arx is None)
        if from_url and (url is not None):
            _id = self.bibcode_from_url(url)
            if _id is not None:
                self._log.info("_id = '{}'".format(_id))
                if self.is_arxiv_id(_id):
                    if arx is None:
                        arx = _id
                        self._log.info("set arx = '{}'".format(arx))
                        self[SOURCE.ARXIVID] = arx
                # if this is not an arxiv ID, it should be a bibcode
                elif bib is None:
                    bib = _id
                    self._log.info("set bib = '{}'".format(bib))
                    self[SOURCE.BIBCODE] = bib

        # Use arxiv to get DOI and url and name (via authors and year)
        #    determine if arxiv is useful here
        from_arx = (doi is None) or (url is None) or (name is None)
        if from_arx and (arx is not None):
            papers = arxiv.query(arx, prune=True, max_results=10)
            if len(papers) == 1:
                papers = papers[0]
                if authors is None:
                    # Select only last names
                    #    e.g. ['C. Uhlemann', 'C. Park', "B. L'Huillier"]
                    authors = [aa.split()[-1] for aa in papers.authors]
                    self._log.info("authors = '{}'".format(authors))
                if year is None:
                    year = papers.published_parsed.tm_year
                    self._log.info("year = '{}'".format(year))
                if url is None:
                    url = papers.id
                    self._log.info("set url = '{}'".format(url))
                    self[SOURCE.URL] = url
                if doi is None:
                    doi = papers.doi
                    self._log.info("set doi = '{}'".format(doi))
                    self[SOURCE.DOI] = doi
            else:
                self._log.warning("Multiple results from arxiv for '{}'".format(arx))

        # use DOI or bib to get the other and name (via authors and year)
        #    determine if doi is useful here
        need_name = ((name is None) and ((authors is None) or (year is None)))
        from_doi = (bib is None) or need_name
        from_bib = (doi is None) or need_name
        if (from_doi and (doi is not None)) or (from_bib and (bib is not None)):
            kw = {'fl': ['bibcode', 'author', 'doi', 'year']}
            if bib is None:
                kw['doi'] = doi
            else:
                kw['bibcode'] = bib
            papers = list(self.ads.SearchQuery(**kw))
            if len(papers) == 1:
                papers = papers[0]
                if bib is None:
                    bib = papers.bibcode
                    self._log.info("set bib = '{}'".format(bib))
                    self[SOURCE.BIBCODE] = bib
                if doi is None:
                    doi = papers.doi
                    self._log.info("set doi = '{}'".format(doi))
                    self[SOURCE.DOI] = doi
                if authors is None:
                    # Select only last names
                    #    e.g. ['Kelley, Luke Zoltan', 'Ramirez-Ruiz, Enrico']
                    authors = [aa.split(',')[0] for aa in papers.author]
                    self._log.info("authors = '{}'".format(authors))
                if year is None:
                    year = papers.year
                    self._log.info("year = '{}'".format(year))
            else:
                self._log.warning("Multiple results from ads for doi: '{}'".format(doi))

        # Construct URL
        if url is None:
            if bib is not None:
                url = "http://adsabs.harvard.edu/abs/" + bib
            elif arx is not None:
                url = "https://arxiv.org/abs/" + arx

            if url is not None:
                self[SOURCE.URL] = url

        # Construct name
        if name is None:
            name = self.name_from_authors_and_year(authors, year)
            self[SOURCE.NAME] = name

        self._log.info("    name = '{}', bib = '{}', arx = '{}', url = '{}', doi = '{}'".format(
            name, bib, arx, url, doi))

        return

    def sort_func(self, key):
        if key == self._KEYS.NAME:
            return 'aaa'
        if key == self._KEYS.BIBCODE:
            return 'aab'
        if key == self._KEYS.ARXIVID:
            return 'aac'
        if key == self._KEYS.ALIAS:
            return 'zzz'
        return key

    def append_sources_from(self, other):
        """`CatDict.append_sources_from` should never be called in `Source`.
        """
        raise RuntimeError("`Source.append_sources_from` called.")

    def is_duplicate_of(self, other):
        """Check if this Source is a duplicate of another.

        Unlike the function in the super class, this method will return True
        if *either* name or bibcode is the same.
        """
        # If these are not the same type, return False
        if type(other) is not type(self):
            return False

        # Go over all expected parameters and check equality of each
        for key in self._KEYS.compare_vals():
            # If only one object has this parameter, not the same
            # This is commented out for sources because two sources are
            # considered the same if they share a name *or* a bibcode
            # if (key in self) != (key in other):
            #     continue
            # If self doesnt have this parameter (and thus neither does), skip
            if key not in self or key not in other:
                continue

            # Now, both objects have the same parameter, compare them
            if self[key] == other[key]:
                return True

        return False

    @classmethod
    def bibcode_from_url(cls, url):
        """Given a URL, try to find the ADS bibcode.

        Currently: only ads or arxiv URLs will work, e.g.

        Returns
        -------
        code : str or 'None'
            The Bibcode if found, otherwise 'None'

        """
        try:
            code = url.split('/abs/')
            code = code[1].strip()
            return code
        except:
            return None

    @classmethod
    def is_arxiv_id(cls, arx):
        """Check the given ID str to see if it matches the pattern for an arxiv ID number.
        """
        if arx is None:
            return False
        match = re.search('^[0-9]{4}\.[0-9]{4,5}', arx)
        retval = True if (match is not None) else False
        return retval

    @classmethod
    def name_from_authors_and_year(cls, authors, year):
        """Create a citation formatted name based on an author list and year.

        Arguments
        ---------
        authors : (N,) str
            Authors should be a list of str, where each entry is the author last-name only.
        year : str or int
            Year of publication, if a string then this should be a four-digit year.

        Returns
        -------
        name : str

        """
        if (authors is None) or (year is None):
            return None

        if len(authors) == 1:
            name = "{} {}".format(authors[0], year)
        elif len(authors) < cls._MAX_AUTHORS_TO_LIST:
            name = "{} & {} {}".format(", ".join(authors[:-1]), authors[-1], year)
        else:
            name = "{} et al. {}".format(authors[0], year)

        return name
