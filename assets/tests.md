# Test different executions of the pipeline

## Normal Run
- regular run: bam, rg --> 20260525
- regular run: split --> 20260525
- regular run: with fixed references --> 20260525
- regular run: with multiple fixed references --> 20260525

## Other TaxLvls
- regular run: with taxlvl order --> 20260525
- regular run: with taxlvl order and (multiple) fixed references --> 20260525
- regular run: with taxlvl genus --> 20260525
- regular run: with taxlvl genus and (multiple) fixed references --> 20260525

## Reruns
- rerun: with fixed families after regular run --> 20260525
- rerun: after multiple fixed references, with multiple fixed reference --> 20260525
- rerun: after taxlvl order with fixed references (default) --> 20260525
- rerun: after taxlvl order with fixed references (and taxlvl o) --> 20260525
- rerun: after taxlvl genus