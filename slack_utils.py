import os

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

    if settings is None:
        return(None)

    getSlackToken(token_file=token_file)
    client = slackclient(os.environ.get('SLACK_TOKEN'))

    return(client)

def postSlackMessage(client, event, channel='slack_bot_test'):

    summary = '<!channel> \n'+event['slack']

    # Post summary data
    kwargs = {}
    kwargs['text'] = summary
    kwargs['channel'] = channel
    kwargs['username'] = settings.botname

    resp = slackCommand(client, 'chat.postMessage', **kwargs)
