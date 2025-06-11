# set_debug_level RSZ rebuffer 3
# set_debug_level RSZ repair_setup 2
puts "reading lib.."
foreach libFile [glob "../ASAP7/LIB/*nldm*.lib"] {
    puts "lib: $libFile"
    read_liberty $libFile
}
puts "reading lef.."
read_lef ../ASAP7/techlef/asap7_tech_1x_201209.lef
foreach lef [glob "../ASAP7/LEF/*.lef"] {
    read_lef $lef
}
puts "reading def.."
read_def ../aes_cipher_top/aes_cipher_top.def
read_sdc ../aes_cipher_top/aes_cipher_top.sdc
source ../ASAP7/setRC.tcl

# write_bookshelf ../aes_cipher_top/aes_cipher_top_rsz_before.pl


estimate_parasitics -placement
report_tns
report_power
report_hpwl

set db [::ord::get_db]
set block [[$db getChip] getBlock]
set insts [$block getInsts]

set inst_positions [dict create]

foreach inst $insts {
    lassign [$inst getLocation] x y
    set inst_name [$inst getName]
    dict set inst_positions $inst_name [list $x $y]
    # puts "inst: $inst_name, ($x, $y)"
}


improve_placement

estimate_parasitics -placement
report_tns
report_power
report_hpwl

repair_timing -skip_pin_swap -skip_gate_cloning -skip_buffer_removal -verbose

estimate_parasitics -placement
report_tns
report_power
report_hpwl

detailed_placement

for {set i 0} {$i < 10} {incr i} {
    improve_placement
}

estimate_parasitics -placement
report_tns
report_power
report_hpwl

puts "\n=== Comparing initial and final positions ==="

set total_displacement 0
set moved_instances 0
set total_instances 0

foreach inst $insts {
    set inst_name [$inst getName]
    lassign [$inst getLocation] final_x final_y
    lassign [dict get $inst_positions $inst_name] initial_x initial_y
    
    set dx [expr abs($final_x - $initial_x)]
    set dy [expr abs($final_y - $initial_y)]
    set distance [expr $dx + $dy]
    
    incr total_instances
    set total_displacement [expr $total_displacement + $distance]
    
    if {$distance > 0} {
        incr moved_instances
        # puts "inst: $inst_name moved from ($initial_x, $initial_y) to ($final_x, $final_y), distance: $distance"
    }
}

puts "\n=== Displacement Statistics ==="
puts "Total instances: $total_instances"
puts "Moved instances: $moved_instances"
puts "Unmoved instances: [expr $total_instances - $moved_instances]"
puts "Total displacement: $total_displacement"

if {$total_instances > 0} {
    set avg_displacement [expr double($total_displacement) / ($total_instances * 1000)]
    puts "Average displacement (all instances): [format "%.2f" $avg_displacement]"
}

write_def ../aes_cipher_top/aes_cipher_top_rsz_after.def
# write_bookshelf ../aes_cipher_top/aes_cipher_top_rsz_after.pl
exit
