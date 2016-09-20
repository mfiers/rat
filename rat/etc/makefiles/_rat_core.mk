
main: all
rmain: rat_refresh all

rat_refresh:
	@echo "Refreshing: $(RAT_TEMPLATE)"
	rat mf install -f $(RAT_TEMPLATE)

