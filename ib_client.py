'''
IBClient class for accessing Interactive Broker data.

@author:  neilswainston
'''
import sys

import Queue as queue
import ib.ext.Contract
import ib.opt


class IBClient(object):
    '''Client class for accessing Interactive Broker data.'''

    def __init__(self, timeout=20, **kwargs):
        self.__queue = queue.Queue()
        self.__timeout = timeout

        self.__connection = ib.opt.ibConnection(**kwargs)
        self.__connection.connect()

    def get_portfolio(self):
        '''Get portfolio data.'''
        self.__connection.reqAccountUpdates(True, 'D999999')

        portfolio = []

        # Listen for responses from self.__connection:
        while True:
            _, _, msg = self.__queue.get(timeout=self.timeout)

            if isinstance(msg, ib.opt.message.accountDownloadEnd):
                # If account data downloaded fully, stop listening:
                break

            if isinstance(msg, ib.opt.message.updatePortfolio):
                contract = msg.contract

                ticker = '%s-%s-%s' % (contract.m_symbol,
                                       contract.m_secType,
                                       contract.m_exchange)

                portfolio.append(ticker)

        return portfolio


def main(args):
    '''Test method.'''
    client_id = args[0]

    try:
        ibm = IBClient(client_id=client_id)
        portfolio = ibm.get_portfolio()
        print portfolio
    except queue.Empty:
        print 'Timeout occured'


if __name__ == '__main__':
    main(sys.argv[1:])
