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
/^class Hjelmslev/ { start_write(0) }
/^\s*def init/ { start_write(0) }
/^\s*def init_spanning_points/ { start_write(0) }
/^\s*def span/ { start_write(0) }
/^\s*def splits/ { start_write(0) }
/^\s*def check_root_group/ { start_write(0) }
/^\s*def partial_root_group_elts/ { start_write(0) }
/^\s*def enough_roots/ { start_write(0) }
/^\s*def chambers_of_some_apartment/ { start_write(0) }
/^\s*def moufang/ { start_write(0) }
/^\s*def gap_graph/ { start_write(0) }
/^\s*def stabilizers/ {start_write(0) }
/^\s*Adjacencies/ { start_write(0) }
{ if (write == 1)
  {
    gsub("    ", "  ");
    print
  }
  count = count - 1
  if (count == 0)
    stop_write()
}
!/^\s*$/ { put_dots() }

