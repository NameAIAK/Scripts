import time
import os
import logging
from datetime import datetime
import sys

## Set up the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)    
logger = logging.getLogger(__name__)
cmd = 'wget -c https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz --no-check-certificate'
while True:
    try:
        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        logger.info(f"{time_current}")
        os.system(cmd) 
        time.sleep(300)
    except:
        time.sleep(300)