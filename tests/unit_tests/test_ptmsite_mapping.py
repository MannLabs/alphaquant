import re
import pytest
import alphaquant.ptm.ptmsite_mapping as aqptm
import pandas as pd

# Constants
INPUT_FILE = "../../test_data/unit_tests/ptmsite_mapping/shortened_aq.tsv.zip"
SAMPLEMAP_FILE = "../../test_data/unit_tests/ptmsite_mapping/samplemap.tsv"

# Helper functions
def compare_parsed_and_handpicked_site_and_prot(refprot_parsed, refprot_handpicked, site_parsed, sites_handpicked):
    if refprot_parsed == refprot_handpicked:
        return compare_sites(site_parsed, sites_handpicked)
    return False

def compare_sites(sites_parsed, sites_handpicked):
    if 'nan' in sites_parsed:
        return False  # we ignore the nan cases in this test
    
    numbers1 = extract_numbers(sites_parsed)
    numbers2 = extract_numbers(sites_handpicked)
    
    numbers1.sort()
    numbers2.sort()

    return numbers1 == numbers2

def extract_numbers(lst):
    if isinstance(lst, str):
        lst = eval(lst)
    return [int(re.search(r'\d+', str(item)).group()) for item in lst]

def get_site_handpicked(protseq, pepseq, ptm_positions_peptide):
    start_idx_protein = protseq.index(pepseq)
    sites_handpicked = str([x+start_idx_protein for x in ptm_positions_peptide])
    return sites_handpicked

# Fixture for the mapped DataFrame
@pytest.fixture(scope="module")
def mapped_df():
    input_df = pd.read_csv(INPUT_FILE, sep="\t", nrows=1000)
    samplemap_df = pd.read_csv(SAMPLEMAP_FILE, sep="\t")
    df = aqptm.assign_dataset(input_df, results_folder=None, samplemap_df=samplemap_df)[0]
    df['site'] = [str(x) for x in df['site']]
    return df

# Test function
def check_df_is_as_expected(mapped_df, refprot_handpicked, sites_handpicked, num_matches_expected):
    num_matches = 0
    for reprot_parsed, site_parsed in zip(mapped_df["REFPROT"], mapped_df["site"]):
        site_parsed = site_parsed.replace("S", "").replace("T", "").replace("Y", "")
        matches = compare_parsed_and_handpicked_site_and_prot(reprot_parsed, refprot_handpicked, site_parsed, sites_handpicked)
        if matches:
            num_matches += 1
    assert num_matches == num_matches_expected, f"Expected {num_matches_expected} matches, but got {num_matches} matches"

# Test cases
def test_p49023_case1(mapped_df):
    refprot = "P49023"
    protseq = "MDDLDALLADLESTTSHISKRPVFLSEETPYSYPTGNHTYQEIAVPPPVPPPPSSEALNGTILDPLDQWQPSSSRFIHQQPQSSSPVYGSSAKTSSVSNPQDSVGSPCSRVGEEEHVYSFPNKQKSAEPSPTVMSTSLGSNLSELDRLLLELNAVQHNPPGFPADEANSSPPLPGALSPLYGVPETNSPLGGKAGPLTKEKPKRNGGRGLEDVRPSVESLLDELESSVPSPVPAITVNQGEMSSPQRVTSTQQQTRISASSATRELDELMASLSDFKIQGLEQRADGERCWAAGWPRDGGRSSPGGQDEGGFMAQGKTGSSSPPGGPPKPGSQLDSMLGSLQSDLNKLGVATVAKGVCGACKKPIAGQVVTAMGKTWHPEHFVCTHCQEEIGSRNFFERDGQPYCEKDYHNLFSPRCYYCNGPILDKVVTALDRTWHPEHFFCAQCGAFFGPEGFHEKDGKAYCRKDYFDMFAPKCGGCARAILENYISALNTLWHPECFVCRECFTPFVNGSFFEHDGQPYCEVHYHERRGSLCSGCQKPITGRCITAMAKKFHPEHFVCAFCLKQLNKGTFKEQNDKPYCQNCFLKLFC"
    pepseq = "FIHQQPQSSSPVYGSSAK"
    ptm_positions_peptide = [10]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=4)

def test_p49023_case2(mapped_df):
    refprot = "P49023"
    protseq = "MDDLDALLADLESTTSHISKRPVFLSEETPYSYPTGNHTYQEIAVPPPVPPPPSSEALNGTILDPLDQWQPSSSRFIHQQPQSSSPVYGSSAKTSSVSNPQDSVGSPCSRVGEEEHVYSFPNKQKSAEPSPTVMSTSLGSNLSELDRLLLELNAVQHNPPGFPADEANSSPPLPGALSPLYGVPETNSPLGGKAGPLTKEKPKRNGGRGLEDVRPSVESLLDELESSVPSPVPAITVNQGEMSSPQRVTSTQQQTRISASSATRELDELMASLSDFKIQGLEQRADGERCWAAGWPRDGGRSSPGGQDEGGFMAQGKTGSSSPPGGPPKPGSQLDSMLGSLQSDLNKLGVATVAKGVCGACKKPIAGQVVTAMGKTWHPEHFVCTHCQEEIGSRNFFERDGQPYCEKDYHNLFSPRCYYCNGPILDKVVTALDRTWHPEHFFCAQCGAFFGPEGFHEKDGKAYCRKDYFDMFAPKCGGCARAILENYISALNTLWHPECFVCRECFTPFVNGSFFEHDGQPYCEVHYHERRGSLCSGCQKPITGRCITAMAKKFHPEHFVCAFCLKQLNKGTFKEQNDKPYCQNCFLKLFC"
    pepseq = "FIHQQPQSSSPVYGSSAK"
    ptm_positions_peptide = [13]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=6)

def test_q15185_case1(mapped_df):
    refprot = "Q15185"
    protseq = "MQPASAKWYDRRDYVFIEFCVEDSKDVNVNFEKSKLTFSCLGGSDNFKHLNEIDLFHCIDPNDSKHKRTDRSILCCLRKGESGQSWPRLTKERAKLNWLSVDFNNWKDWEDDSDEDMSNFDRFSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEKMPDLE"
    pepseq = "FSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEK"
    ptm_positions_peptide = [29]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=4)

def test_q15185_case2(mapped_df):
    refprot = "Q15185"
    protseq = "MQPASAKWYDRRDYVFIEFCVEDSKDVNVNFEKSKLTFSCLGGSDNFKHLNEIDLFHCIDPNDSKHKRTDRSILCCLRKGESGQSWPRLTKERAKLNWLSVDFNNWKDWEDDSDEDMSNFDRFSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEKMPDLE"
    pepseq = "FSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEK"
    ptm_positions_peptide = [26]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=1)

def test_q15185_case3(mapped_df):
    refprot = "Q15185"
    protseq = "MQPASAKWYDRRDYVFIEFCVEDSKDVNVNFEKSKLTFSCLGGSDNFKHLNEIDLFHCIDPNDSKHKRTDRSILCCLRKGESGQSWPRLTKERAKLNWLSVDFNNWKDWEDDSDEDMSNFDRFSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEKMPDLE"
    pepseq = "FSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEK"
    ptm_positions_peptide = [26, 29]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=11)

def test_q15185_case4(mapped_df):
    refprot = "Q15185"
    protseq = "MQPASAKWYDRRDYVFIEFCVEDSKDVNVNFEKSKLTFSCLGGSDNFKHLNEIDLFHCIDPNDSKHKRTDRSILCCLRKGESGQSWPRLTKERAKLNWLSVDFNNWKDWEDDSDEDMSNFDRFSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEKMPDLE"
    pepseq = "FSEMMNNMGGDEDVDLPEVDGADDDSQDSDDEK"
    ptm_positions_peptide = [2, 26, 29]
    sites_handpicked = get_site_handpicked(protseq, pepseq, ptm_positions_peptide)
    check_df_is_as_expected(mapped_df, refprot_handpicked=refprot, sites_handpicked=sites_handpicked, num_matches_expected=2)