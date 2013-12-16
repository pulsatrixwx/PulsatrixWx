#!/usr/bin/perl

#
#	parse_stations.pl
#	Purpose:	Extract station location information from a remote data source
#	Started:	19 June 2010 by Tim Supinie (tsupinie@ou.edu)
# 	Completed:	19 June 2010
#	Modified:	[not yet]
#

# Remote data file that will be used as the source
my $remote_data_file_name = "http://www.rap.ucar.edu/weather/surface/stations.txt";

# Input file name (file name from $remote_data_file_name)
my $file_name = "stations.txt";

# Surface station output file name
my $sfc_file_name = "sfc_stations.csv";

# Upper air station output file name
my $ua_file_name = "ua_stations.csv";

# Array of character positions to break each line at
my @breaks = (0, 3, 20, 26, 32, 39, 42, 47, 50, 55, 62, 65, 68, 71, 74, 77, 79, 81, 83);

# Work array
my @work =   ("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "");

# Holding arrays for surface and upper air stations
my @sfc_stations = ();
my @ua_stations = ();

# Retrieve remote data file
system "wget $remote_data_file_name > /dev/null 2>&1";

# Read input file
open(STATIONFILE, "<", $file_name) or die "Can't open $file_name: $!\n";
while (<STATIONFILE>) {
	# Read each line
	if (substr($_, 0, 1) ne "!" and length() > 80) {
		# The line is not a comment and its length is more than 80 characters
		#   (filter out things that aren't possible reporting points).
		for (my $idx = 0; $idx < $#breaks; $idx++) {
			# Parse out value
			$work[$idx] = substr $_, $breaks[$idx], $breaks[$idx + 1] - $breaks[$idx];

			# Strip leading and trailing white space
			$work[$idx] =~ s/^[\s]*([\S]*(?:[\s][\S]+)*)[\s]*$/$1/mg
		}

		# Convert coordinates in degrees/minutes/seconds to decimal degrees
		$min_latitude = substr $work[6], 0, length($work[6]) - 1;
		$hem_latitude = substr $work[6], -1, 1;
		$latitude = ($work[5] + $min_latitude / 60.0) * ($hem_latitude eq "N" ? 1 : -1);

		$min_longitude = substr $work[8], 0, length($work[8]) - 1;
		$hem_longitude = substr $work[8], -1, 1;
		$longitude = ($work[7] + $min_longitude / 60.0) * ($hem_longitude eq "E" ? 1 : -1);

		if (($work[17] eq "US" or $work[17] eq "CA" or $work[17] eq "MX") and $work[10] eq "X") {
			# Keep US, Canadian, and Mexican METAR reporting stations
			push @sfc_stations, "$work[0],$work[1],$work[2],$latitude,$longitude,$work[9]";
		}

		if (($work[17] eq "US" or $work[17] eq "CA" or $work[17] eq "MX") and $work[13] eq "X") {
			# Keep US, Canadian, and Mexican upper air stations
			push @ua_stations, "$work[0],$work[1],$work[2],$latitude,$longitude,$work[9]";
		}
	}
}

close STATIONFILE;

# Write out surface stations
open(SFCFILE, ">", $sfc_file_name) or die "Can't open $sfc_file_name: $!\n";
foreach (@sfc_stations) {
	print SFCFILE "$_\n";
}
close SFCFILE;

# Write out upper air stations
open(UAFILE, ">", $ua_file_name) or die "Can't open $ua_file_name: $!\n";
foreach (@ua_stations) {
	print UAFILE "$_\n";
}
close UAFILE;

# Clean up (delete input file)
unlink $file_name;
