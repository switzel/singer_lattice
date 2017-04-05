function stop_write()
{
  if (write == 1) {
    write = 0
    dots = 0
  }
}

function start_write(n)
{
  if (write == 0)
    if (dots == 1)
      print "[...]"
    else
      print ""
  write = 1
  count = n
}

function put_dots()
{
  dots = 1
}

BEGIN { stop_write() }
END { start_write() }
/^$/ { stop_write() }
/^diff_sets/ { start_write(5) }
/^\s*def eq_classes/ { start_write(0) }
/^class DiffSetFactory/ { start_write(0) }
/^\s*def difference_matrices/ { start_write(0) }
/^\s*def aut_diff_set/ { start_write(4) }
/^\s*def transform/ { start_write(0) }
/^\s*def diff_mat_diff_set/ {start_write(0) }
/^\s*def orbit/ { start_write(0) }
/^\s*dsf = DiffSetFactory/ { start_write(4) }
/^\s*if __name__/ { start_write(1) }
{ if (write == 1)
    print
  count = count - 1
  if (count == 0)
    stop_write()
}
!/^\s*$/ { put_dots() }

