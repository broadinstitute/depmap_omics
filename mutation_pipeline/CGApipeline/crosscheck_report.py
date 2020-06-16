import csv
import argparse

def main():

	parser = argparse.ArgumentParser(description='Parse the CrosscheckFingerprints metrics file produced by GATK4.0.4.0 and make a report.')
	parser.add_argument("metrics_file", type=argparse.FileType('r'), help="CrosscheckFingerprints output metrics file.")
	args = parser.parse_args()


	metrics = csv.DictReader(args.metrics_file, delimiter="\t")
	
	min_lod = 1000
	min_lod_left_lane = ""
	min_lod_right_lane = ""

	metrics_dict = {}

	#Go over the metrics file and populate metrics_dict with the LOD results per pair of lanes
	for met in metrics:
		met_left_lane = met["LEFT_FILE"].split("/")[-1] + ":" + met["LEFT_GROUP_VALUE"]
		met_right_lane = met["RIGHT_FILE"].split("/")[-1] + ":" + met["RIGHT_GROUP_VALUE"]
		met_lod = float(met["LOD_SCORE"])

		if met_lod < min_lod:
			min_lod = met_lod
			min_lod_left_lane = met_left_lane
			min_lod_right_lane = met_right_lane

		if met_left_lane not in metrics_dict:
			metrics_dict[met_left_lane] = {}
			metrics_dict[met_left_lane][met_right_lane] = met_lod

		else:
			metrics_dict[met_left_lane][met_right_lane] = met_lod

	#Write the minimum LOD value and the two lanes associated with it in to .txt files to be read by FireCloud later.
	with open("crosscheck_min_lod_value.txt", "w") as f:
		f.write(str(min_lod))

	with open("crosscheck_min_lod_lanes.txt", "w") as f:
		f.write("{0}, {1}".format(min_lod_left_lane, min_lod_right_lane))


	#Start creating the report file
	report_file = open("report.html", "w")

	#Write the used style definitions into the html file
	report_file.write("<style>\n")
	report_file.write("\t<!--\n")
	report_file.write("\t\ttd.good {\n")
	report_file.write("\t\t\tbackground-color: #B2F7AB\n")
	report_file.write("\t\t}\n")
	report_file.write("\t\ttd.warning {\n")
	report_file.write("\t\t\tbackground-color: #FCCB60\n")
	report_file.write("\t\t}\n")
	report_file.write("\t\ttd.bad {\n")
	report_file.write("\t\t\tbackground-color: #F78691\n")
	report_file.write("\t\t}\n")
	report_file.write("\t\ttd.self {\n")
	report_file.write("\t\t\tbackground-color: #C9C7C8\n")
	report_file.write("\t\t}\n")
	report_file.write("\t-->\n")
	report_file.write("</style>\n\n")

	#Write the file header
	report_file.write("<h3>Lane Crosscheck QC</h3>\n\n")

	#Write the overall result of the crosschecking process
	if min_lod < 1:
		report_file.write("<br><font color=\"red\">PROBLEMS DETECTED</font><br>\n\n")
	elif min_lod < 15:
		report_file.write("<br><font color=\"orange\">SUSPICIOUSLY LOW CONCORDANCE</font><br>\n\n")
	else:
		report_file.write("<br><font color=\"green\">ALL LANES OK</font><br>\n\n")
	report_file.write("<hr>\n\n")
	
	#Create a sorted list of all lanes to be used in order to create the LOD score table.
	all_lanes = sorted(list(metrics_dict.keys()))

	#Start the table by writing a row of all column headers
	report_file.write("<table>\n\t<tr>\n\t\t<th></th>\n")
	for lane in all_lanes:
		report_file.write("\t\t<th class=\"vertical\">{0}</th>\n".format(lane))

	report_file.write("\t</tr>\n\n")

	#Populate the table
	report_row_counter = 0
	for lane in all_lanes:
		report_file.write("\t<tr>\n")
		report_file.write("\t\t<th>{0}</th>\n".format(lane))

		for i in xrange(report_row_counter):
			column_lane = all_lanes[i]
			try:
				lod_score = metrics_dict[lane][column_lane]
			except:
				lod_score=-999
			if lod_score >= 15:
				result_class = "good"
			elif lod_score >= 1:
				result_class = "warning"
			else:
				result_class = "bad"
			report_file.write("\t\t<td class=\"{0}\">{1}</td>\n".format(result_class, lod_score))

		report_file.write("\t\t<td class=\"self\"></td>\n")

		remaining_rows = len(all_lanes) - report_row_counter - 1
		for i in xrange(remaining_rows):
			report_file.write("\t\t<td></td>\n")

		report_file.write("\t</tr>\n")
		report_row_counter += 1

	report_file.write("</table>\n")


if __name__ == '__main__':
	main()








