import logging
from contextlib import contextmanager

LOG_TEMPLATE_ROW, LOG_TEMPLATE_CHAIN = None, None
LOG_FILE_SEPARATOR = " |LOG_FILE_DETAILS: "

class MemoryHandler(logging.Handler):
    log_records = []
    errors_per_row = {"CHAIN1": {}, "CHAIN2": {}}

    def emit(self, record):
        self.log_records.append({
            "row": str(LOG_TEMPLATE_ROW),
            "chain": str(LOG_TEMPLATE_CHAIN),
            "time": self.formatter.formatTime(record),
            "level": record.levelname,
            "msg": record.getMessage()
        })

        self._add_row_error(record)

    def _add_row_error(self, record):
        if LOG_TEMPLATE_ROW is not None:
            mssg = f"CHAIN {LOG_TEMPLATE_CHAIN} {record.levelname}: " + record.getMessage()
            if LOG_FILE_SEPARATOR in mssg:
                mssg = mssg.split(LOG_FILE_SEPARATOR)[0]

            if int(LOG_TEMPLATE_ROW) not in self.errors_per_row[f"CHAIN{LOG_TEMPLATE_CHAIN}"]:
                self.errors_per_row[f"CHAIN{LOG_TEMPLATE_CHAIN}"][int(LOG_TEMPLATE_ROW)] = mssg
            else:
                self.errors_per_row[f"CHAIN{LOG_TEMPLATE_CHAIN}"][int(LOG_TEMPLATE_ROW)] += " | "
                self.errors_per_row[f"CHAIN{LOG_TEMPLATE_CHAIN}"][int(LOG_TEMPLATE_ROW)] += mssg


@contextmanager
def log_with_error_collector():
    logger = logging.getLogger()
    memory_handler = MemoryHandler()
    logger.addHandler(memory_handler)

    try:
        yield memory_handler
    finally:
        logger.removeHandler(memory_handler)
        memory_handler.close()


def format_log_with_context(logging_info):
    for file_handler in logging.root.handlers:
        file_handler.setFormatter(logging.Formatter(f"%(asctime)s %(levelname)s {logging_info}: %(message)s"))

def format_log_neutral():
    for file_handler in logging.root.handlers:
        file_handler.setFormatter(logging.Formatter(f"%(asctime)s %(levelname)s: %(message)s"))

@contextmanager
def set_logging_context(template_row=None, chain=None, log_context=None):
    global LOG_TEMPLATE_ROW
    global LOG_TEMPLATE_CHAIN

    LOG_TEMPLATE_ROW = str(template_row) if template_row is not None else None
    LOG_TEMPLATE_CHAIN = str(chain) if chain is not None else None

    if log_context is None:
        if template_row is not None:
            if chain is not None:
                log_context = f"ROW {template_row} CHAIN {chain}"
            else:
                log_context = f"ROW {template_row}"
        else:
            assert False, "template_row or log_context must be set"

    format_log_with_context(log_context)
    try:
        yield
    finally:
        LOG_TEMPLATE_ROW = None
        LOG_TEMPLATE_CHAIN = None

        format_log_neutral()
