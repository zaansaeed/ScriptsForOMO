from abnumber import Chain

seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'
chain = Chain(seq, scheme='imgt')

# chain.print(numbering=True)
# print(chain.cdr3_seq)
# chain.print(numbering=True)
# chain.print_tall()

for pos, aa in chain:
    print(pos, aa)