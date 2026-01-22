import os
import requests

def slackCommand(client, command, **kwargs):
    if command in ['conversations.list']:
        resp = client.api_call(command, params=kwargs)
    else:
        resp = client.api_call(command, json=kwargs)
    message = '#' * 80 + '\n'
    message += '#' * 80 + '\n'
    message += f'################ SLACK COMMAND {command} ################\n\n\n'
    print(message)
    return(resp)

def getSlackToken(token_file=None):
    if token_file and os.path.exists(token_file):
        with open(token_file) as f:
            token = f.readline().replace('\n', '')
            os.environ['SLACK_TOKEN'] = token
    else:
        print(f'WARNING: token file {token_file} does not exist.')
        print('Create token file or manually set SLACK_TOKEN variable.')

def setUpSlack(token_file):
    from slack_sdk import WebClient as slackclient
    getSlackToken(token_file=token_file)
    client = slackclient(os.environ.get('SLACK_TOKEN'))
    return(client)

def postSlackMessage(client, message, channel='slack_bot_test',
    botname='Robotron', at_channel=True):
    
    if at_channel:
        summary = '<!channel> \n'+message
    else:
        summary = message

    # Post summary data
    kwargs = {}
    kwargs['text'] = summary
    kwargs['channel'] = channel
    kwargs['username'] = botname
    resp = slackCommand(client, 'chat.postMessage', **kwargs)

def getFileUploadURL(client, file_path, channel):
    file_name = os.path.basename(file_path)
    file_size = os.path.getsize(file_path)

    response_get_url = client.files_getUploadURLExternal(
            filename=file_name,
            length=file_size
        )
    upload_url = response_get_url["upload_url"]
    file_id = response_get_url["file_id"]

    with open(file_path, "rb") as f:
        file_content = f.read()
        
    headers = {"Content-Type": "application/octet-stream"}
    response_upload = requests.post(upload_url, data=file_content, headers=headers)
    response_upload.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

    return(file_id)


def postSlackFiles(client, file_ids, message, channel, at_channel=True):

    if at_channel:
        summary = '<!channel> \n'+message
    else:
        summary = message

    resp = slackCommand(client, 'conversations.list')
    channel_id = [r['id'] for r in resp['channels'] if r['name_normalized']==channel]

    if len(channel_id)==0:
        return(None)
    else:
        CHANNEL_ID = channel_id[0]

    response_complete_upload = client.files_completeUploadExternal(
        files=[{"id": f} for f in file_ids],
        channel_id=CHANNEL_ID,
        initial_comment=summary
    )

    if response_complete_upload["ok"]:
        print(f"File '{file_ids}' uploaded and shared successfully!")
    else:
        print(f"Error completing upload: {response_complete_upload['error']}")
