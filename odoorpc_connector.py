#!/usr/bin/env python3
# Author(s): Silvio Gregorini (silviogregorini@openforce.it)
# Copyright 2019 Silvio Gregorini (silviogregorini@openforce.it)
# License LGPL-3.0 or later (https://www.gnu.org/licenses/lgpl).

from odoorpc import ODOO


class Odoo(ODOO):

    def __init__(
            self,
            db, login, pw,
            host='localhost', protocol='jsonrpc',
            port=8069, timeout=120,
            version=None, opener=None
    ):
        # For privacy and security reasons, we will not hardcode values
        # needed for connecting to any database.
        # Whoever runs this script, will have to know the right inputs to
        # make it work properly.
        if not all([db, login, pw]):
            raise ValueError("Please specify database, login and password")
        super().__init__(host, protocol, port, timeout, version, opener)
        self.login(db, login, pw)

    def get_env(self):
        return self._env

    def get_model(self, class_name):
        env = self.get_env()
        return env[class_name] if class_name in env else ''

    def get_record(self, class_name, res_ids):
        if isinstance(res_ids, int):
            res_ids = [res_ids]
        return self.get_model(class_name).browse(res_ids)
