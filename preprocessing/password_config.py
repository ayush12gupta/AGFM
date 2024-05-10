import json, os


DIR_PATH = os.path.realpath(os.path.dirname(__file__))
with open(f'{DIR_PATH}/configs/credentials.json', 'r') as f:
        config_cred = json.load(f)


unavuser=config_cred["UNAV_user"]
unavpass=config_cred["UNAV_user"]

asfuser=config_cred["ASF_user"]
asfpass=config_cred["ASF_password"]

eossouser=config_cred["EOSSO_user"]
eossopass=config_cred["EOSSO_user"]

