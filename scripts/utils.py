"""
This file contains some functions that are used by several other scripts.
"""

ALLOWED_CHARS = "ACDEFGHIKLMNPQRSTVWY[]1459."

def find_modifications(peptide, style):
	"""Receives a peptide sequence where the modifications are formatted like: KIM[15.9949]AKEMI
			Returns the peptide and its modifications either in the style of pDeep or Prosit.
			Note that currently, only oxidized Methionine is supported as a PTM.
			Also, only standard amino acids are supported, so e.g. all peptides with Selenocystein are dropped"""
	if any(char not in ALLOWED_CHARS for char in peptide):
		return None
	if style.lower() == "pdeep":
		res = ""
		pos = peptide.find("M[")
		while pos != -1:
			res += "{},Oxidation[M];".format(pos + 1)
			peptide = peptide[:peptide.find("[")] + peptide[peptide.find("]") + 1:]
			pos = peptide.find("M[")
		assert "[" not in peptide, "Illegal Modification in Peptide: {}".format(peptide)
		return peptide, res
	elif style.lower() == "prosit":
		peptide = peptide.replace("M[15.9949]", "M(ox)")
		assert "[" not in peptide, "Illegal Modification in Peptide: {}".format(peptide)
		return peptide


def convert_seq_notation(sequence, source="pdeep", target="prosit"):
	"""Convertes between different notation styles of peptide sequences with modifications. Currently only supports
		conversion from pdeep to prosit notation."""
	if source == "pdeep" and target == "prosit":
		seq, modifications = sequence.split("|")[:2]
		offset = 0
		for mod in modifications.split(";")[:-1]:
			pos, mod_name = mod.split(",")
			pos = int(pos)
			assert mod_name == "Oxidation[M]", "pdeep sequence title contains illegal modification"
			assert seq[pos + offset - 1] == 'M', "pdeep spectrum title contains modification at non-Methionine residue"
			seq = seq[:pos + offset] + "(ox)" + seq[pos + offset:]
			offset += 4  # Every insertion shifts the sequence by 4 (len('(ox)'))
		return seq
	else:
		return sequence


def add_rt_to_spectra(msms_mgf, rt_df, index_list):
	res = []
	for i in index_list:
		spectrum = msms_mgf.get_by_id(i)
		spectrum['params']['irt'] = rt_df[
			rt_df.modified_sequence == convert_seq_notation(
				spectrum['params']['title'])].iRT.values[0]
		res.append(spectrum)
	return res