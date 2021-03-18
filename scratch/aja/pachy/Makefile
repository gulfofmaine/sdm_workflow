# Docker bits
up: down
	docker-compose up -d --build
	docker-compose logs -f
	
down:
	docker-compose down

logs:
	docker-compose logs -f

prune:
	docker volume rm $(shell docker volume ls -qf dangling=true)
	docker system prune -a

#####
## Workflow pipeline steps
# Sample data prep
samp_dat_prepped.rds: 
	./pipelines/samp_dat_prep/R/samp_dat_prep.R "./pfs/raw_dat/" 1 1985 2005 "./pfs/samp_dat_prep/"
