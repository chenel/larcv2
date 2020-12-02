#!/usr/bin/env python3
import ROOT, sys
from ROOT import std, TChain
from larcv import larcv

if len(sys.argv) < 2:
   print ('Usage: python',sys.argv[0],'CONFIG_FILE [LARCV_FILE1 LARCV_FILE2 ...]')
   sys.exit(1)

# Instantiate, configure, and ensure it's kWrite mode
proc = larcv.ProcessDriver('ProcessDriver')
proc.configure(sys.argv[1])
if not proc.io().io_mode() == 1:
	print('Supera needs to be run with IOManager IO mode kWrite! exiting...')
	sys.exit(1)

# Set up input
ch = TChain('EDepSimEvents')
if len(sys.argv) > 2:
	for argv in sys.argv[2:]:
		if not argv.endswith('.root'): continue
		print('Adding input:',argv)
		ch.AddFile(argv)
print("Chain has", ch.GetEntries(), "entries")
event_range = (proc.batch_start_entry(), proc.batch_start_entry() + proc.batch_num_entry()) \
	if proc.batch_num_entry() > 0 \
	else (0, ch.GetEntries())
print('Processing', event_range[1] - event_range[0], 'events')
sys.stdout.flush()

# Initialize and retrieve a list of processes that belongs to SuperaBase inherited module classes
proc.initialize()
supera_procs = []
for name in proc.process_names():
	pid = proc.process_id(name)
	module = proc.process_ptr(pid)
	if getattr(module,'is')('Supera'):
		print('Running a Supera module:',name)
		supera_procs.append(pid)

# Event loop
for entry in range(*event_range):
	print("considering event:", entry)
	sys.stdout.flush()
	bytes = ch.GetEntry(entry)
	if bytes < 1:
		break

	ev = ch.Event 

	proc.set_id(ev.RunId,0,ev.EventId)

	# set event pointers
	for pid in supera_procs:
		module = proc.process_ptr(pid)
		module.SetEvent(ev)

	proc.process_entry()
proc.finalize()