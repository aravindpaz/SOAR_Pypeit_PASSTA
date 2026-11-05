"""Slack Web API helpers for uploading spectra and formatted messages."""

import logging
import os
from typing import Any

logger = logging.getLogger(__name__)

import requests


def slack_command(client: Any, command: str, **kwargs: Any) -> Any:
    """Invoke a Slack Web API method and print a debug banner.

    Parameters
    ----------
    client : slack_sdk.web.client.WebClient
        Authenticated Slack client.
    command : str
        API method name (e.g. ``'chat.postMessage'``).
    **kwargs
        Forwarded as JSON body or query parameters depending on the method.

    Returns
    -------
    Any
        Parsed API response from ``client.api_call``.
    """
    if command in {'conversations.list'}:
        resp = client.api_call(command, params=kwargs)
    else:
        resp = client.api_call(command, json=kwargs)
    message = '#' * 80 + '\n'
    message += '#' * 80 + '\n'
    message += f'################ SLACK COMMAND {command} ################\n\n\n'
    logger.debug('%s', message.rstrip())
    return resp


def get_slack_token(token_file: str | None = None) -> None:
    """Read a bot token from disk and set ``SLACK_TOKEN`` in the environment.

    Parameters
    ----------
    token_file : str or None
        Path to a one-line file containing the token. If missing, prints a
        warning and leaves the environment unchanged.
    """
    if token_file and os.path.exists(token_file):
        with open(token_file) as f:
            token = f.readline().replace('\n', '')
            os.environ['SLACK_TOKEN'] = token
    else:
        logger.warning('Token file %s does not exist.', token_file)
        logger.warning('Create token file or manually set SLACK_TOKEN variable.')


def set_up_slack(token_file: str) -> Any:
    """Create a Slack ``WebClient`` using a token file.

    Parameters
    ----------
    token_file : str
        Path passed to :func:`get_slack_token`.

    Returns
    -------
    slack_sdk.web.client.WebClient
        Configured client using ``SLACK_TOKEN`` from the environment.
    """
    from slack_sdk import WebClient as slackclient

    get_slack_token(token_file=token_file)
    return slackclient(os.environ.get('SLACK_TOKEN'))


def post_slack_message(
    client: Any,
    message: str,
    channel: str = 'slack_bot_test',
    botname: str = 'Robotron',
    at_channel: bool = True,
) -> Any:
    """Post a plain-text message to a channel.

    Parameters
    ----------
    client : slack_sdk.web.client.WebClient
        Authenticated Slack client.
    message : str
        Message body.
    channel : str, optional
        Channel name (default ``slack_bot_test``).
    botname : str, optional
        Legacy display username for the post.
    at_channel : bool, optional
        If True, prefix with ``<!channel>``.

    Returns
    -------
    Any
        Slack API response from ``chat.postMessage``.
    """
    summary = ('<!channel> \n' + message) if at_channel else message
    kwargs = {'text': summary, 'channel': channel, 'username': botname}
    return slack_command(client, 'chat.postMessage', **kwargs)


def get_file_upload_url(client: Any, file_path: str, channel: str) -> str:
    """Stage a file with Slack's external upload flow and return its file id.

    Parameters
    ----------
    client : slack_sdk.web.client.WebClient
        Authenticated Slack client.
    file_path : str
        Local path to the file to upload.
    channel : str
        Unused in this step; reserved for future channel scoping.

    Returns
    -------
    str
        Slack ``file_id`` for use with :func:`post_slack_files`.

    Notes
    -----
    ``channel`` is kept for API compatibility with older call sites.
    """
    _ = channel
    file_name = os.path.basename(file_path)
    file_size = os.path.getsize(file_path)

    response_get_url = client.files_getUploadURLExternal(
        filename=file_name,
        length=file_size,
    )
    upload_url = response_get_url['upload_url']
    file_id = response_get_url['file_id']

    with open(file_path, 'rb') as f:
        file_content = f.read()

    headers = {'Content-Type': 'application/octet-stream'}
    response_upload = requests.post(upload_url, data=file_content, headers=headers)
    response_upload.raise_for_status()

    return file_id


def post_slack_files(
    client: Any,
    file_ids: list[str],
    message: str,
    channel: str,
    at_channel: bool = True,
) -> None:
    """Complete an external upload and share files to a channel.

    Parameters
    ----------
    client : slack_sdk.web.client.WebClient
        Authenticated Slack client.
    file_ids : list of str
        File ids from :func:`get_file_upload_url`.
    message : str
        Initial comment shown with the upload.
    channel : str
        Human-readable channel name (resolved via ``conversations.list``).
    at_channel : bool, optional
        If True, prefix the comment with ``<!channel>``.
    """
    summary = ('<!channel> \n' + message) if at_channel else message

    resp = slack_command(client, 'conversations.list')
    channel_id = [r['id'] for r in resp['channels'] if r['name_normalized'] == channel]

    if len(channel_id) == 0:
        return

    channel_slack_id = channel_id[0]

    response_complete_upload = client.files_completeUploadExternal(
        files=[{'id': f} for f in file_ids],
        channel_id=channel_slack_id,
        initial_comment=summary,
    )

    if response_complete_upload['ok']:
        logger.info('Slack upload completed for file ids: %s', file_ids)
    else:
        logger.error(
            'Slack upload failed: %s',
            response_complete_upload.get('error', 'unknown'),
        )
