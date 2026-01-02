#!/usr/bin/env sh
################################################################################
#
# pomerol2triqs
#
# Copyright (C) 2017-2026 Igor Krivenko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

#
# Update copyright years in files
#

FIND=$(which find)
SED=$(which sed)

# Collect paths of files to be updated
FILES="README.md"

for PATTERN in "*.cpp" "*.hpp" "*.py" "*.in" "CMakeLists.txt"
do
    FILES="$FILES $($FIND . -name $PATTERN)"
done

# Write the new year string
YEARS_REGEX="([0-9]{4})-[0-9]{4}"
YEARS_UPDATED="\\1-$(date +%Y)"

COPYRIGHT_REGEX="Copyright \(C\) ${YEARS_REGEX}"
COPYRIGHT_UPDATED="Copyright \(C\) ${YEARS_UPDATED}"
for FILE in $FILES
do
    $SED -i -E "s/${COPYRIGHT_REGEX}/${COPYRIGHT_UPDATED}/g" $FILE
done
