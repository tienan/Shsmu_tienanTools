#haow

awk -F, '{print $1","$4","$7}' ukb9820.csv | sed -e 's/3-2.0/sex/g' -e  's/eid/id/g' | head 